# Preliminary analysis script for CH in HCT project
# Loads and filters initial somatic variant calls, assigns quality metrics to samples
# Last edited: February 2024

#------------------------------------------------------------------------------
# Required Libraries
#------------------------------------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(readxl)         # Reading Excel files
library(janitor)        # Data cleaning functions

#------------------------------------------------------------------------------
# External Data Referenced (not explicitly loaded)
#------------------------------------------------------------------------------
# - /Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/MergedResults/merge_cD3.maf
#   Contains raw somatic variant calls from initial pipeline
#
# - Data_in/hsmetrics.csv
#   Contains sequencing QC metrics from MultiQC output
#
# - Per-sample start/stop duplicate rate files at:
#   /Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/cD3/{sample}/results/QC/{sample}.SSduprate.hist.txt

#------------------------------------------------------------------------------
# Load and Process Metadata
#------------------------------------------------------------------------------
# Clinical sample metadata
SequencingMetadataMastersheet <- read_excel(
  "~/Library/CloudStorage/OneDrive-MontefioreMedicine/CH in HCT/database/SequencingMetadataMastersheet.xlsx", 
  na = c("NA","")) %>% 
  filter(Exclude != 1) %>% 
  mutate_at("Group", as.factor)

# Extract key metadata for analysis
MetadataIDs <- SequencingMetadataMastersheet %>%
  select(Library_ID, Group, Sample_ID, SampleType, SampleTimepoint, CellPopulation, isScreen) %>%
  mutate(Time = case_when(
    SampleTimepoint == "Donor Day 1" ~ 0,
    SampleTimepoint == "Donor Pre"  ~ -10,
    SampleTimepoint == "one year" ~ 365,
    SampleTimepoint == "Day 150" ~ 150,
    SampleTimepoint == "Day 100" ~ 100,
    SampleTimepoint == "Day 75" ~ 75,
    SampleTimepoint == "Day 180" ~ 180,
    SampleTimepoint == "Day 50" ~ 50,
    SampleTimepoint == "Day 75-90 Pre 1st DLI"  ~ 75
  ))

# Clinical data and demographics
DemoDat <- read_excel("~/Library/CloudStorage/OneDrive-MontefioreMedicine/CH in HCT/database/Data for Jacob 3.2024 - De-IDed.xlsx", 
                      col_types = c("text", "numeric", "numeric", "text", "skip", 
                                  "skip", "numeric", "text", "numeric", "text", 
                                  "date", "text", "text", "text", "text", 
                                  "text", "text", "text", "numeric", "date", 
                                  "text", "text", "text", "text", "numeric", 
                                  "date", "skip", "skip")) %>% 
  clean_names(case = "upper_camel", abbreviations = c("ID", "CH", "CD")) %>% 
  filter(!is.na(DonorSampleID))

#------------------------------------------------------------------------------
# Load and Process Variant Data
#------------------------------------------------------------------------------
# Define variant classifications considered non-synonymous
nonsyn <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
            "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", 
            "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")

# Load raw somatic variant calls
RawSomatic_MAF <- read_delim("/Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024//MergedResults/merge_cD3.maf",
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
  mutate(ConsensusDepth = "cD3") %>%
  select(where(~!all(is.na(.x)))) %>%
  mutate(VAF = t_alt_count / t_depth,
         Variant_Class = if_else(Variant_Classification %in% nonsyn, 
                               "non_synonymous", "synonymous"))

#------------------------------------------------------------------------------
# Process Library Quality Metrics
#------------------------------------------------------------------------------
# Count passing variants per library
PASSperLib <- RawSomatic_MAF %>% 
  filter(FILTER == "PASS") %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  summarise(n = n())

# Load and process hybrid selection metrics
hsmetrics <- read_csv("Data_in/hsmetrics.csv") %>%
  mutate(Library_ID = str_split_i(`Sample Name`, pattern = " | ", i = 3))

# Process start/stop duplicate rates
DupRates <- lapply(SequencingMetadataMastersheet$Library_ID,
                  function(x) read_table(paste0("/Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/cD3/",
                                              x,"/results/QC/",x,".SSduprate.hist.txt"), 
                                       col_names = c(x,"DupRate"))) %>% 
  reduce(full_join, by = "DupRate")

# Summarize duplicate rates and combine with other metrics
DupRatesSummary <- DupRates %>% 
  pivot_longer(cols = starts_with("CiH")) %>% 
  group_by(name) %>% 
  summarise(total = sum(value, na.rm = TRUE),
            single = value[DupRate == 1],
            double = value[DupRate == 2]) %>% 
  mutate(rate = 1-single/total,
         dupover2 = 1-(single+double)/total) %>% 
  left_join(SequencingMetadataMastersheet %>% select(Library_ID, LibInput), 
            by = c("name" = "Library_ID")) %>% 
  left_join(hsmetrics %>% select(Library_ID, `Mean target coverage`), 
            by = c("name" = "Library_ID")) %>% 
  left_join(PASSperLib, by = c("name" = "Tumor_Sample_Barcode"))

#------------------------------------------------------------------------------
# Identify Low Quality Libraries
#------------------------------------------------------------------------------
# Flag libraries failing quality thresholds
LQlibs <- DupRatesSummary %>% 
  filter(LibInput < 25 | dupover2 > 0.05) %>% 
  pull(name)

# Identify groups with low quality screen samples
LQgroups <- SequencingMetadataMastersheet %>% 
  filter(Library_ID %in% LQlibs & isScreen == 1) %>%
  pull(Group) %>% unique()

#------------------------------------------------------------------------------
# Filter Somatic Variants
#------------------------------------------------------------------------------
# Get passing variant positions from raw calls
vcf_pos.pass <- RawSomatic_MAF$vcf_pos[RawSomatic_MAF$FILTER == "PASS"]

# Apply comprehensive filtering
MAF.pass <- RawSomatic_MAF %>% 
  filter(vcf_pos %in% vcf_pos.pass) %>%
  left_join(SequencingMetadataMastersheet, by = c("Tumor_Sample_Barcode" = "Library_ID")) %>% 
  mutate(VarID = paste0(Group,Hugo_Symbol,HGVSc,Chromosome,Start_Position)) %>%
  group_by(vcf_pos) %>%
  mutate(filt.MultipleGroups2 = length(unique(Group)) > 2) %>% 
  filter(Tumor_Sample_Barcode %in% MetadataIDs$Library_ID) %>%
  filter(!is.na(Group)) %>%
  group_by(VarID) %>%
  mutate(
    # Technical filters
    filt.LowAlt = max(t_alt_count) < 5,
    filt.LowCoverage = max(t_depth) < 600,
    filt.LowVAF = max(VAF) < 0.0025,
    filt.SingleLQSample = length(unique(Tumor_Sample_Barcode)) == 1 & 
      Tumor_Sample_Barcode %in% LQlibs,
    # Biological filters
    filt.Germline = grepl("rs", dbSNP_RS) & any(str_detect(FILTER,"germline")),
    filt.AllFail = !any(str_detect(FILTER,"PASS")),
    # Final pass/fail assignment
    PASS = !if_any(.cols = starts_with("filt.")))

# Write filtered MAF for review
MAF.pass %>% write_csv("Data_explore/MAF.pass.csv")

#------------------------------------------------------------------------------
# Identify Validated Variants
#------------------------------------------------------------------------------
# Get variants validated in screen samples
VARsIDedSCREEN <- MAF.pass %>% 
  filter(VAF > 0.005,
         PASS == TRUE,
         FILTER == "PASS",
         isScreen == 1) %>%
  pull(VarID) %>% 
  unique()

# Add manually validated variant
VARsIDedSCREEN <- c(VARsIDedSCREEN, "AS062909TET2c.3765C>Gchr4105243740")

# Get variants validated in donor samples
VARsIDedDONOR <- MAF.pass %>% 
  filter(VAF > 0.005,
         PASS == TRUE,
         FILTER == "PASS",
         SampleType == "DONOR") %>%
  pull(VarID) %>% 
  unique()

# Save key data for other scripts
save(SequencingMetadataMastersheet, MetadataIDs, DemoDat, LQlibs, LQgroups,
     VARsIDedSCREEN, VARsIDedDONOR, GroupsWPass,
     file = "RData/01.needed.RData")

save(RawSomatic_MAF, DupRatesSummary, hsmetrics, MAF.pass,
     file = "RData/01.other.RData")