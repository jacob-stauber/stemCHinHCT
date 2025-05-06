# Force Calling Script for CH in HCT project
# Prepares input files for force calling variants across sample groups
# Variants initially detected in any sample are force-called across all samples in a group
# Last edited: February 2024

#------------------------------------------------------------------------------
# Required Libraries
#------------------------------------------------------------------------------
# Required libraries are inherited from 01.Prelim.R:
# - tidyverse: Data manipulation and visualization
# - readxl: Reading Excel files

#------------------------------------------------------------------------------
# External Data Referenced (not explicitly loaded)
#------------------------------------------------------------------------------
# From 01.Prelim.R:
# - RawSomatic_MAF
# - DupRatesSummary
# - hsmetrics
# - MAF.pass
# - SequencingMetadataMastersheet
# - MetadataIDs
# - DemoDat
# - LQlibs
# - LQgroups
# - VARsIDedSCREEN
# - VARsIDedDONOR

#------------------------------------------------------------------------------
# Global Variables
#------------------------------------------------------------------------------
# Base path for force calling results
FCpath = "/Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/ForceCalls/01.indivFC/"

#------------------------------------------------------------------------------
# Prepare Sample Groups for Force Calling
#------------------------------------------------------------------------------
# Get groups with passing variants
GroupsWPass <- MAF.pass %>% 
  filter(PASS == TRUE) %>% 
  pull(Group) %>% 
  unique()

# Create BAM file list for force calling
# One line per sample: Group,BAM_path
SequencingMetadataMastersheet %>% 
  filter(Group %in% GroupsWPass) %>% 
  select(Group, Library_ID) %>% 
  mutate(Library_ID = paste0(
    "/gs/gsfs0/home/jstaube1/CHinHCT/RESULTS/FINAL/B_APRIL2024/cD3/",
    Library_ID,
    "/results/bam/",
    Library_ID,
    ".bam"
  )) %>%
  write_csv(paste0(FCpath, "filelist_cD3.indiv2.txt"), 
            col_names = FALSE)

#------------------------------------------------------------------------------
# Functions for MAF File Generation
#------------------------------------------------------------------------------
#' Create MAF Files for Force Calling
#'
#' Creates filtered MAF files per group for use in force calling.
#' NOTE: Function is study-specific, check variables if rerunning.
#'
#' @param GroupID Group identifier
#' @param MafFilt Filtered MAF data frame (default: MAF.pass)
#' @param path Output directory path
#' @return None, writes MAF file to disk
MAFsForForceCalls <- function(GroupID, 
                             MafFilt = MAF.pass, 
                             path = paste0(FCpath,"inMAFs/")) {
  MafFilt %>% 
    filter(Group == GroupID) %>%
    filter(PASS == TRUE) %>%
    ungroup() %>%
    distinct(VarID, .keep_all = TRUE) %>%
    select(any_of(colnames(RawSomatic_MAF))) %>%
    write_delim(paste0(path, GroupID, ".sites.maf"), 
                delim = "\t")
}

#------------------------------------------------------------------------------
# Generate Force Calling Input Files
#------------------------------------------------------------------------------
# Create MAF files for each group
sapply(GroupsWPass, MAFsForForceCalls)

#------------------------------------------------------------------------------
# Save Data for Downstream Analysis
#------------------------------------------------------------------------------
# Save core data needed for other scripts
save(SequencingMetadataMastersheet, MetadataIDs, DemoDat, LQlibs, LQgroups,
     VARsIDedSCREEN, VARsIDedDONOR, GroupsWPass,
     file = "RData/01.needed.RData")

# Save additional analysis data
save(RawSomatic_MAF, DupRatesSummary, hsmetrics, MAF.pass,
     file = "RData/01.other.RData")

#------------------------------------------------------------------------------
# Next Steps
#------------------------------------------------------------------------------
# Run ForceCall.sh on HPC cluster for variant calling
# Key steps performed in shell script:
# 1. Convert MAF to VCF format for Mutect2 input:
#    conda activate vcf2maf
#    maf2vcf.pl \
#      --input-maf ./inMAFs/{GROUP}.sites.maf \
#      --output-dir {GROUP} \
#      --ref-fasta $genome_file
#
# 2. Sort and index VCF:
#    conda activate TSA
#    vcf-sort {GROUP}/{GROUP}.sites.vcf > {GROUP}/{GROUP}.sites.sort.vcf
#    bgzip {GROUP}/{GROUP}.sites.sort.vcf
#    tabix {GROUP}/{GROUP}.sites.sort.vcf.gz
#
# 3. Force call variants for each sample in group:
#    gatk Mutect2 \
#      -R $genome_file \
#      -L {GROUP}/{GROUP}.sites.sort.vcf.gz \
#      -alleles {GROUP}/{GROUP}.sites.sort.vcf.gz \
#      -O {GROUP}/{GROUP}.{SAMPLE}.FC.vcf \
#      -I {SAMPLE}.bam \
#      --germline-resource $germline_resource
#
# 4. Split multi-allelic sites:
#    gatk LeftAlignAndTrimVariants \
#      -R $genome_file \
#      -V {GROUP}/{GROUP}.{SAMPLE}.FC.vcf \
#      -O {GROUP}/{GROUP}.{SAMPLE}.FC.split.vcf \
#      --split-multi-allelics
#
# 5. Convert back to MAF format:
#    conda activate vcf2maf
#    vcf2maf.pl \
#      --input-vcf {GROUP}/{GROUP}.{SAMPLE}.FC.split.vcf \
#      --output-maf {GROUP}/{GROUP}.{SAMPLE}.FC.split.MAF \
#      --ref-fasta $genome_file \
#      --cache-version 105 \
#      --tumor-id {SAMPLE}
#
# 6. Merge group results:
#    head -2 {GROUP}/{GROUP}.{SAMPLE1}.FC.split.MAF > {GROUP}/{GROUP}.FC.cD3.MAF
#    tail -n +3 -q {GROUP}/*.FC.split.MAF >> {GROUP}/{GROUP}.FC.cD3.MAF
#
# Results will be used in 03.FinalVars.R for final variant filtering