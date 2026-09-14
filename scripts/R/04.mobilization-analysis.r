#==============================================================================
# 04.mobilization-analysis.r
#
# Donor pre- vs post-mobilization clonal dynamics.
#
# Computes purity-corrected CD34+ variant allele frequencies, propagates the
# correction variance, and classifies each variant as expanding, shrinking, or
# unchanged in the CD34+ and CD34- compartments by two-proportion Z-test.
#
# Calculation only: figure generation is not included.
#
# Inputs:
#   - Final non-synonymous variant table (donor pre- and post-mobilization rows)
#   - Library metadata with CD34 enrichment (flow purity) per library
#
# Outputs (written to output_dir):
#   - Final_ZTest_Significance_Table.csv     full statistics, all columns
#   - Raw_Detection_Dynamics_0.5pct.csv      uncorrected detection calls
#   - Supplementary_Table_Clonal_Dynamics.csv  formatted table for publication
#==============================================================================

library(tidyverse)

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
variants_file  <- "data/FinalVars_nonsynonymous.csv"
metadata_file  <- "data/MetadataIDs.csv"
output_dir     <- "results"

# Timepoint codes for the paired donor samples
time_pre       <- -10   # pre-mobilization draw
time_post      <- 0     # day 1 apheresis, post G-CSF

# Analysis thresholds
vaf_floor           <- 0.001   # floor for corrected VAF (avoids zero variance)
vaf_ceiling         <- 0.500   # biological max for a heterozygous somatic variant
min_max_vaf         <- 0.005   # variant must reach this in at least one fraction
detection_threshold <- 0.005   # threshold for the uncorrected detection calls
p_threshold         <- 0.10    # significance cutoff for the Z-test
default_purity_pct  <- 90      # assumed CD34 purity when flow data is missing

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#------------------------------------------------------------------------------
# Data ingestion
#------------------------------------------------------------------------------
variants_post    <- read_csv(variants_file, show_col_types = FALSE)
library_metadata <- read_csv(metadata_file, show_col_types = FALSE)

#------------------------------------------------------------------------------
# Build the paired pre/post variant table
#------------------------------------------------------------------------------

# Donors with both a pre- and post-mobilization sample in both fractions
valid_donors <- library_metadata %>%
  filter(Time %in% c(time_pre, time_post)) %>%
  filter(CellPopulation %in% c("CD34pos", "CD34neg")) %>%
  group_by(Group) %>%
  summarize(
    has_pre  = any(Time == time_pre),
    has_post = any(Time == time_post),
    .groups  = "drop"
  ) %>%
  filter(has_pre & has_post) %>%
  pull(Group)

mobilization_variants <- variants_post %>%
  filter(Group %in% valid_donors) %>%
  filter(Time %in% c(time_pre, time_post)) %>%
  filter(CellPopulation %in% c("CD34pos", "CD34neg")) %>%
  left_join(
    library_metadata %>% select(Library_ID, CD34Enrichment),
    by = "Library_ID"
  ) %>%
  mutate(CD34Enrichment = ifelse(is.na(CD34Enrichment), default_purity_pct, CD34Enrichment)) %>%
  # Require solid read support at a minimum of one timepoint
  group_by(Group, VarID) %>%
  filter(any(t_alt_count > 1, na.rm = TRUE)) %>%
  ungroup()

# One row per variant, with VAF, depth, alt count and purity for each
# fraction/timepoint combination
master_table <- mobilization_variants %>%
  select(Group, VarID, Hugo_Symbol, HGVSp_Short, Time, CellPopulation,
         VAF, t_depth, t_alt_count, CD34Enrichment) %>%
  pivot_wider(
    id_cols     = c(Group, VarID, Hugo_Symbol, HGVSp_Short),
    names_from  = c(CellPopulation, Time),
    values_from = c(VAF, t_depth, t_alt_count, CD34Enrichment),
    names_sep   = "_"
  ) %>%
  rename(
    Raw_Pre  = !!paste0("VAF_CD34pos_", time_pre),
    Raw_Post = !!paste0("VAF_CD34pos_", time_post),
    Neg_Pre  = !!paste0("VAF_CD34neg_", time_pre),
    Neg_Post = !!paste0("VAF_CD34neg_", time_post),

    D_pos_pre  = !!paste0("t_depth_CD34pos_", time_pre),
    D_pos_post = !!paste0("t_depth_CD34pos_", time_post),
    D_neg_pre  = !!paste0("t_depth_CD34neg_", time_pre),
    D_neg_post = !!paste0("t_depth_CD34neg_", time_post),

    Alt_pos_pre  = !!paste0("t_alt_count_CD34pos_", time_pre),
    Alt_pos_post = !!paste0("t_alt_count_CD34pos_", time_post),
    Alt_neg_pre  = !!paste0("t_alt_count_CD34neg_", time_pre),
    Alt_neg_post = !!paste0("t_alt_count_CD34neg_", time_post),

    Purity_Pre_Pct  = !!paste0("CD34Enrichment_CD34pos_", time_pre),
    Purity_Post_Pct = !!paste0("CD34Enrichment_CD34pos_", time_post)
  ) %>%
  # Missing or non-positive VAFs are set to the floor so downstream variance
  # and log-ratio terms remain defined
  mutate(across(c(Raw_Pre, Raw_Post, Neg_Pre, Neg_Post),
                ~ ifelse(is.na(.) | . <= 0, vaf_floor, .))) %>%
  mutate(Max_VAF = pmax(Raw_Pre, Raw_Post, Neg_Pre, Neg_Post, na.rm = TRUE)) %>%
  filter(Max_VAF >= min_max_vaf)

#------------------------------------------------------------------------------
# Purity correction
#
# The immunomagnetically selected CD34+ fraction contains a residual mature
# cell population. The observed CD34+ VAF is therefore a mixture of the true
# HSPC VAF and the CD34- VAF, weighted by flow-determined purity:
#
#   VAF_corrected = [VAF_CD34pos - (1 - purity) * VAF_CD34neg] / purity
#------------------------------------------------------------------------------
corrected_table <- master_table %>%
  mutate(
    f_pre  = Purity_Pre_Pct / 100,
    f_post = Purity_Post_Pct / 100,

    raw_calc_corr_pre  = (Raw_Pre  - ((1 - f_pre)  * Neg_Pre))  / f_pre,
    raw_calc_corr_post = (Raw_Post - ((1 - f_post) * Neg_Post)) / f_post,

    # Bound the corrected values: over-subtraction can drive the estimate
    # negative, and no heterozygous somatic clone exceeds 50%
    Corr_Pre = case_when(
      raw_calc_corr_pre < vaf_floor   ~ vaf_floor,
      raw_calc_corr_pre > vaf_ceiling ~ vaf_ceiling,
      TRUE                            ~ raw_calc_corr_pre
    ),
    Corr_Post = case_when(
      raw_calc_corr_post < vaf_floor   ~ vaf_floor,
      raw_calc_corr_post > vaf_ceiling ~ vaf_ceiling,
      TRUE                             ~ raw_calc_corr_post
    )
  ) %>%
  select(-raw_calc_corr_pre, -raw_calc_corr_post)

#------------------------------------------------------------------------------
# Variance propagation and Z-test classification
#
# CD34- variants are tested on the observed VAFs directly. CD34+ variants are
# tested on the corrected VAFs, with the contamination term's variance carried
# through the correction formula by the delta method.
#------------------------------------------------------------------------------
significance_table <- corrected_table %>%
  mutate(
    # Binomial sampling variance of each observed VAF
    var_pos_pre  = (Raw_Pre  * (1 - Raw_Pre))  / D_pos_pre,
    var_neg_pre  = (Neg_Pre  * (1 - Neg_Pre))  / D_neg_pre,
    var_pos_post = (Raw_Post * (1 - Raw_Post)) / D_pos_post,
    var_neg_post = (Neg_Post * (1 - Neg_Post)) / D_neg_post,

    # --- CD34- (background) ---
    z_score_neg = (Neg_Post - Neg_Pre) / sqrt(var_neg_post + var_neg_pre + 1e-10),
    p_value_neg = 2 * (1 - pnorm(abs(z_score_neg))),
    log2FC_neg  = log2(Neg_Post / Neg_Pre),
    Class_Neg = case_when(
      p_value_neg >= p_threshold ~ "Unchanged",
      log2FC_neg > 0             ~ "Expanding",
      log2FC_neg < 0             ~ "Shrinking",
      TRUE                       ~ "Unchanged"
    ),

    # --- CD34+ (purity-corrected) ---
    # Delta method on VAF_corr = (VAF_pos - (1-f)*VAF_neg) / f
    var_corr_pre  = ((1 / f_pre)^2  * var_pos_pre)  + (((1 - f_pre)  / f_pre)^2  * var_neg_pre),
    var_corr_post = ((1 / f_post)^2 * var_pos_post) + (((1 - f_post) / f_post)^2 * var_neg_post),

    z_score_pos = (Corr_Post - Corr_Pre) / sqrt(var_corr_post + var_corr_pre + 1e-10),
    p_value_pos = 2 * (1 - pnorm(abs(z_score_pos))),
    log2FC_pos  = log2(Corr_Post / Corr_Pre),
    Class_Pos = case_when(
      p_value_pos >= p_threshold ~ "Unchanged",
      log2FC_pos > 0             ~ "Expanding",
      log2FC_pos < 0             ~ "Shrinking",
      TRUE                       ~ "Unchanged"
    )
  ) %>%
  select(
    # Variant identifiers
    Group, VarID, Hugo_Symbol, HGVSp_Short,

    # Purity and read support
    Purity_Pre_Pct, Purity_Post_Pct, f_pre, f_post,
    D_pos_pre, Alt_pos_pre, D_neg_pre, Alt_neg_pre,
    D_pos_post, Alt_pos_post, D_neg_post, Alt_neg_post,

    # Observed and corrected VAFs
    Raw_Pre, Raw_Post, Neg_Pre, Neg_Post, Corr_Pre, Corr_Post,

    # CD34+ (corrected) results
    log2FC_pos, p_value_pos, Class_Pos, z_score_pos,
    var_corr_pre, var_corr_post, var_pos_pre, var_pos_post,

    # CD34- (background) results
    log2FC_neg, p_value_neg, Class_Neg, z_score_neg,
    var_neg_pre, var_neg_post,

    everything()
  ) %>%
  arrange(p_value_pos)

write_csv(significance_table, file.path(output_dir, "Final_ZTest_Significance_Table.csv"))

#------------------------------------------------------------------------------
# Uncorrected detection dynamics
#
# Companion analysis on the observed VAFs alone, classifying whether each
# variant was detectable above threshold before and after mobilization.
#------------------------------------------------------------------------------
detection_sub_analysis <- significance_table %>%
  mutate(
    pos_pre_det  = Raw_Pre  >= detection_threshold,
    pos_post_det = Raw_Post >= detection_threshold,
    neg_pre_det  = Neg_Pre  >= detection_threshold,
    neg_post_det = Neg_Post >= detection_threshold,

    # Detected in either fraction
    any_pre_det  = pos_pre_det  | neg_pre_det,
    any_post_det = pos_post_det | neg_post_det,

    Pos_Class = case_when(
      pos_pre_det  &  pos_post_det ~ "Present in Both",
      !pos_pre_det &  pos_post_det ~ "Appeared",
      pos_pre_det  & !pos_post_det ~ "Disappeared",
      TRUE                         ~ "Below Threshold"
    ),
    Neg_Class = case_when(
      neg_pre_det  &  neg_post_det ~ "Present in Both",
      !neg_pre_det &  neg_post_det ~ "Appeared",
      neg_pre_det  & !neg_post_det ~ "Disappeared",
      TRUE                         ~ "Below Threshold"
    ),
    Combined_Class = case_when(
      any_pre_det  &  any_post_det ~ "Present in Both",
      !any_pre_det &  any_post_det ~ "Appeared",
      any_pre_det  & !any_post_det ~ "Disappeared",
      TRUE                         ~ "Below Threshold"
    )
  ) %>%
  select(
    Group, VarID, Hugo_Symbol, HGVSp_Short,
    Raw_Pre, Raw_Post, Neg_Pre, Neg_Post,
    Pos_Class, Neg_Class, Combined_Class
  ) %>%
  arrange(Group, Hugo_Symbol)

write_csv(detection_sub_analysis, file.path(output_dir, "Raw_Detection_Dynamics_0.5pct.csv"))

#------------------------------------------------------------------------------
# Formatted supplementary table
#------------------------------------------------------------------------------
supp_table_export <- significance_table %>%
  mutate(
    Variant = paste(Hugo_Symbol, str_remove(HGVSp_Short, "p."), sep = " ")
  ) %>%
  select(
    Group,
    Variant,

    `CD34neg Pre VAF`  = Neg_Pre,
    `CD34neg Post VAF` = Neg_Post,
    `CD34neg Log2FC`   = log2FC_neg,
    `CD34neg P-Value`  = p_value_neg,

    `CD34pos Pre VAF (Corrected)`  = Corr_Pre,
    `CD34pos Post VAF (Corrected)` = Corr_Post,
    `CD34pos Log2FC`               = log2FC_pos,
    `CD34pos P-Value`              = p_value_pos,

    `Purity Pre (%)`  = Purity_Pre_Pct,
    `Purity Post (%)` = Purity_Post_Pct
  ) %>%
  mutate(
    across(c(`CD34neg Pre VAF`, `CD34neg Post VAF`,
             `CD34pos Pre VAF (Corrected)`, `CD34pos Post VAF (Corrected)`),
           ~ round(., 4)),
    across(c(`CD34neg Log2FC`, `CD34pos Log2FC`), ~ round(., 2)),
    across(c(`CD34neg P-Value`, `CD34pos P-Value`),
           ~ formatC(., format = "e", digits = 2)),
    across(c(`Purity Pre (%)`, `Purity Post (%)`), ~ round(., 1))
  ) %>%
  arrange(Group, Variant)

write_csv(supp_table_export, file.path(output_dir, "Supplementary_Table_Clonal_Dynamics.csv"))

#==============================================================================
# End
#==============================================================================
