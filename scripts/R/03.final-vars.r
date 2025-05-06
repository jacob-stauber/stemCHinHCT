# Final Variant Processing Script
# Processes force-called variants, remaps group IDs, and adds metadata
# Last edited: February 2025

#------------------------------------------------------------------------------
# Required Libraries
#------------------------------------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(data.table)     # Fast data manipulation
library(janitor)        # Data cleaning functions

#------------------------------------------------------------------------------
# Load Source Data
#------------------------------------------------------------------------------
# Load previous analysis data
load("RData/01.needed.RData")
load("RData/01.other.RData")

# Load chimerism data
ChimerismData <- read_csv("Data_in/ChimerismData.csv")

# Load group ID mapping
group_mapping <- read_csv("Data_in/groupid_mapping.csv")

#------------------------------------------------------------------------------
# Process Final Variants
#------------------------------------------------------------------------------
# Base path for force calling results
FCpath = "/Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/ForceCalls/01.indivFC/"

# Read force called MAFs for each group
FC <- lapply(GroupsWPass, function(Group) {
  read_delim(paste0(FCpath, Group, "/", Group, ".FC.cd3.MAF"),
             delim = "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE, 
             skip = 1)
})

# Combine and process force called variants
FC_MAF <- data.table::rbindlist(FC) %>% 
  # Add basic annotations
  mutate(
    ConsensusDepth = "cD3",
    Variant_Class = if_else(Variant_Classification %in% nonsyn, 
                            "non_synonymous", "synonymous"),
    VAF = t_alt_count / t_depth,
    # Clean allele encoding
    across(c(Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2), 
           ~gsub("TRUE", "T", .x)),
    # Create identifiers
    VarID = paste0(Group, Hugo_Symbol, HGVSc, Chromosome, Start_Position),
    VCFVAR = paste(Chromosome, Start_Position, Reference_Allele, 
                   Tumor_Seq_Allele2, sep = "-"),
    # Add validation flags
    IN_DONOR = VarID %in% VARsIDedDONOR,
    IN_SCREEN = VarID %in% VARsIDedSCREEN
  ) %>%
  left_join(SequencingMetadataMastersheet, 
            by = c("Tumor_Sample_Barcode" = "Library_ID"))

#------------------------------------------------------------------------------
# Add Additional Annotations and Metrics
#------------------------------------------------------------------------------
# Load pathogenic annotations
annoPD <- read_csv("Data_in/ReportedAsCHinFC.anno.0404.csv")

# Process variants with additional metrics
FC_MAF_processed <- FC_MAF %>%
  group_by(VarID) %>%
  mutate(
    # Calculate group-level metrics
    VAFinDONOR = mean(VAF[SampleType == "DONOR"]),
    maxVAF = max(VAF),
    nPASS = sum(OrigFilter == "PASS", na.rm = TRUE),
    nPresent = sum(VAF > 0),
    nCount = sum(t_alt_count > 1),
    nCountDonor = sum(t_alt_count[SampleType == "DONOR"] > 1)
  ) %>%
  # Join annotations
  left_join(annoPD %>% select(VCFVAR, PD), by = "VCFVAR") %>%
  left_join(MAF.pass %>% select(VarID, Tumor_Sample_Barcode, FILTER) %>% 
              rename(OrigFilter = FILTER))


#------------------------------------------------------------------------------
# Add Additional Annotations and Metrics
#------------------------------------------------------------------------------
# Load pathogenic annotations
# CH-PD annotations from curated database of CH variants
# PD: Pathogenic/likely pathogenic designation based on published criteria
annoPD <- read_csv("Data_in/ReportedAsCHinFC.anno.0404.csv")

# Process variants with additional metrics
# Calculate per-variant metrics for review and filtering
FC_MAF_processed <- FC_MAF %>%
  group_by(VarID) %>%
  mutate(
    # Calculate group-level metrics
    VAFinDONOR = mean(VAF[SampleType == "DONOR"]),
    maxVAF = max(VAF),
    nPASS = sum(OrigFilter == "PASS", na.rm = TRUE),
    nPresent = sum(VAF > 0),
    nCount = sum(t_alt_count > 1),
    nCountDonor = sum(t_alt_count[SampleType == "DONOR"] > 1)
  ) %>%
  # Join annotations
  left_join(annoPD %>% select(VCFVAR, PD), by = "VCFVAR") %>%
  left_join(MAF.pass %>% select(VarID, Tumor_Sample_Barcode, FILTER) %>% 
              rename(OrigFilter = FILTER))

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Prepare Variants for Manual Review
#------------------------------------------------------------------------------
# Generate review file for variants requiring IGV validation
FC_MAF_for_review <- FC_MAF_processed %>%
  select(Group, SampleType, SampleTimepoint, CellPopulation, 
         Tumor_Sample_Barcode, VCFVAR, Hugo_Symbol, HGVSc, HGVSp_Short,
         ConsensusDepth, t_depth, t_alt_count, VAF, Variant_Class,
         Variant_Classification, Variant_Type, VAFinDONOR, IN_DONOR, 
         IN_SCREEN, OrigFilter, PD, Existing_variation, LibInput, 
         isScreen, VarID, maxVAF, nPASS, nPresent, nCount, nCountDonor)

# Write review file
write_csv(FC_MAF_for_review, "Data_explore/FC_MAF_for_review.csv")

#------------------------------------------------------------------------------
# Manual Review Process
#------------------------------------------------------------------------------
# Variants are reviewed in IGV and classified according to these criteria:
# Reason definitions:
#   Confirmed: Variants detected in multiple independent libraries or validated by IGV review
#   FC fail: Variants called in initial pipeline but not reproduced in force calling
#   Germline: Allele frequency > 0.4 and not in CH hotspot
#   Haplotype: Failed IGV verification due to complex haplotype from local realignment
#   Haplotype Dup: Confirmed variants in same haplotype as another; more deleterious retained
#   Missing: Mutect2 fails to call locus in most samples in group
#   Outlier Gene: SAMD9 and SAMD9L calls with complex haplotype patterns
#   Outlier Library: Indels unique to libraries with overrepresented INDEL calls
#   Recip Germline: Allele frequency pattern matches recipient germline vs chimerism
#   Single Template: All variant-supporting reads share same start-stop positions
#   Weak Evidence Indel: Indel >5bp with only one PASS and no spanning reads
#   Below Threshold: All samples below 0.005 allele frequency

# Load manually reviewed variants
FC_MAF_reviewed <- read_csv("Data_in/FC_MAF_reviewed.csv")

#------------------------------------------------------------------------------
# Create Final Variant List
#------------------------------------------------------------------------------
#' Structure of final variant data frame (FinalVars):
#' @field Group              Sample group identifier
#' @field isScreen           Whether sample is screening (1) or not (0)
#' @field SampleType         DONOR or RECIPIENT
#' @field SampleTimepoint    Timepoint of sample collection
#' @field Time               Numeric time relative to transplant (days)
#' @field CellPopulation     Cell type of sample
#' @field Library_ID         Library identifier
#' @field VarID              Unique variant identifier
#' @field VCFVAR             Variant in CHR-POS-REF-ALT format
#' @field Hugo_Symbol        Gene name
#' @field HGVSc              Coding change
#' @field HGVSp_Short        Protein change
#' @field t_depth            Total read depth
#' @field t_alt_count        Alternative allele count
#' @field VAF                Variant allele frequency (alternate reads/total reads)
#' @field Variant_Classification Specific variant type
#' @field Variant_Type       SNP, INS, or DEL
#' @field isNonSyn           Synonymous or non_synonymous
#' @field PD                 Putative CH driver designation
#' @field Existing_variation Known variant IDs
#' @field OrigFilter         Original filter status
#' @field IN_DONOR           Present in donor sample
#' @field IN_SCREEN          Identified in screening
#' @field DonorChimerism     Donor chimerism in recip sample
#' @field ChimerismSD        Donor chimerism SD

# Add chimerism and time information
FinalVars <- FC_MAF_reviewed %>%
  filter(grepl("Confirmed", Reason)) %>%
  # Add time information
  mutate(Time = case_when(
    SampleTimepoint == "Donor Day 1" ~ 0,
    SampleTimepoint == "Donor Pre" ~ -10,
    SampleTimepoint == "one year" ~ 365,
    SampleTimepoint == "Day 150" ~ 150,
    SampleTimepoint == "Day 100" ~ 100,
    SampleTimepoint == "Day 75" ~ 75,
    SampleTimepoint == "Day 180" ~ 180,
    SampleTimepoint == "Day 50" ~ 50,
    SampleTimepoint == "Day 75-90 Pre 1st DLI" ~ 75
  )) %>%
  # Add chimerism data
  left_join(ChimerismData, by = c("Library_ID")) %>%
  mutate(DonorChimerism = 1-Chimerism,
         ChimerismSD = SD) %>% 
  # Add new group IDs
  left_join(group_mapping %>% select(OldID, NewID, ScreenType), 
            by = c("Group" = "OldID")) %>%
  rename(OldGroup = Group, Group = NewID) %>%
  # Split into non-synonymous and all variants
  mutate(isNonSyn = Variant_Class == "non_synonymous") %>%
  select(Group, OldGroup, ScreenType,
         isScreen, SampleType, SampleTimepoint, Time, CellPopulation,
         Library_ID, VarID, VCFVAR, Hugo_Symbol, HGVSc, HGVSp_Short,
         t_depth, t_alt_count, VAF, Variant_Classification, Variant_Type, isNonSyn,
         PD, Existing_variation, OrigFilter, IN_DONOR, IN_SCREEN,
         DonorChimerism, ChimerismSD)

# Create split data lists
FinalVars_split <- list(
  nonsynonymous = list(
    POST = FinalVars %>% 
      filter(isNonSyn == TRUE, ScreenType == "POST") %>%
      select(-isNonSyn),
    PRE = FinalVars %>% 
      filter(isNonSyn == TRUE, ScreenType == "PRE") %>%
      select(-isNonSyn)
  ),
  all = list(
    POST = FinalVars %>% 
      filter(ScreenType == "POST") %>%
      select(-isNonSyn),
    PRE = FinalVars %>% 
      filter(ScreenType == "PRE") %>%
      select(-isNonSyn)
  )
)

#------------------------------------------------------------------------------
# Save Results
#------------------------------------------------------------------------------
save(FinalVars_split, file = "RData/final_variants_split.RData")


write_csv(FinalVars_split$nonsynonymous$POST, "Results/FinalVars_nonsynonymous_POST.csv")
write_csv(FinalVars_split$nonsynonymous$PRE, "Results/FinalVars_nonsynonymous_PRE.csv")
write_csv(FinalVars_split$all$POST, "Results/FinalVars_all_POST.csv")
write_csv(FinalVars_split$all$PRE, "Results/FinalVars_all_PRE.csv")
