# Germline Analysis and Chimerism Calculation
# Processes joint called variants to calculate chimerism
# Last edited: February 2024

#------------------------------------------------------------------------------
# Required Libraries
#------------------------------------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(vcfR)          # VCF file handling
library(pheatmap)      # Heatmap visualization
library(broom)         # Statistical output formatting
library(GGally)        # Correlation plots

#------------------------------------------------------------------------------
# External Data Referenced (not explicitly loaded)
#------------------------------------------------------------------------------
# - /Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/germline/CiH_ALL.recalibrated.funcotator.maf
#   Contains joint-called germline variants annotated by funcotator
#
# - MetadataIDs from previous analysis containing sample grouping information

#------------------------------------------------------------------------------
# Load and Process Germline Variants
#------------------------------------------------------------------------------
# Read joint-called variants
germlineMAF.func <- read_delim(
  "/Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/germline/CiH_ALL.recalibrated.funcotator.maf",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE, comment = "##", skip = 1
)

# Clean and process variants
germlineMAF.func.clean <- germlineMAF.func %>%
  # Remove unknown/problematic entries
  filter(Hugo_Symbol != "__UNKNOWN__") %>%
  filter(!grepl("\\|", gnomAD_genome_AF)) %>%
  # Clean columns
  select(where(~!all(is.na(.x)))) %>%
  select(where(~!all(.x == "__UNKNOWN__"))) %>%
  type_convert() %>%
  # Fix DEL positions
  mutate(vcfpos = ifelse(Variant_Type == "DEL", 
                        Start_Position-1, 
                        Start_Position))

# Find common alleles (>2% frequency in any gnomAD population)
CommonAlleles <- germlineMAF.func.clean %>% 
  filter(if_any(.cols = c(gnomAD_genome_AF_popmax,
                         gnomAD_exome_AF_popmax), 
                ~ . > 0.02)) %>%
  select(vcfpos) %>% 
  unlist()

#------------------------------------------------------------------------------
# Process VCF for Common Variants
#------------------------------------------------------------------------------
# Read and process VCF
germlinevcf <- read.vcfR("/Volumes/users/Jacob Stauber/CHinHCT/RESULTS/FINAL/B_APRIL2024/germline/CiH_ALL.recalibrated.dbsnp.vcf")
germlinevcf.tidy <- vcfR2tidy(germlinevcf, format_fields = c("GT", "AD", "DP"))

#' Calculate VAF from Allelic Depth
#' @param AD Allelic depth string "ref,alt"
#' @return Numeric VAF
VAFfromAD <- function(AD) {
  str_split_fixed(AD, ",", 2) %>%
    apply(1, function(x) as.numeric(x[2])/sum(as.numeric(x)))
}

# Process common alleles
CommonAllelesGT <- germlinevcf.tidy$fix[,1:8] %>%
  full_join(germlinevcf.tidy$gt) %>%
  filter(FILTER == "PASS") %>%
  filter(POS %in% CommonAlleles) %>%
  filter(!grepl(",", ALT)) %>%
  mutate(VAF = VAFfromAD(gt_AD)) %>%
  left_join(MetadataIDs, by = c("Indiv" = "Library_ID"))

#------------------------------------------------------------------------------
# Sample Swap Check Functions
#------------------------------------------------------------------------------
# Prepare annotation for plots
annotation_col <- MetadataIDs %>%
  select(Library_ID, Group, Sample_ID, CellPopulation) %>%
  column_to_rownames(var = "Library_ID")

#' Check for Sample Swaps within a Group
#' @param GroupSelect Group identifier
#' @param AlleleTable Table of common alleles (default: CommonAllelesGT)
#' @param annotation Sample annotations (default: annotation_col)
#' @param exclude Samples to exclude
#' @return Generates heatmap plot
SampleSwapCheck <- function(GroupSelect, 
                          AlleleTable = CommonAllelesGT,
                          annotation = annotation_col,
                          exclude = NA) {
  mat <- AlleleTable %>%
    filter(!Indiv %in% exclude) %>%
    filter(Group == GroupSelect) %>%
    select(ID, Indiv, VAF) %>%
    filter(!is.na(ID)) %>% 
    pivot_wider(names_from = Indiv, values_from = VAF) %>%
    column_to_rownames(var = "ID") %>%
    as.matrix() %>%
    na.omit()
  
  mat <- mat[rowSums(mat) > 0,]
  
  pheatmap(mat, 
           annotation_col = annotation,
           show_rownames = FALSE,
           filename = paste0("Figs/sampleswap/",
                           GroupSelect,
                           ".sampleswapcheck.png"),
           cellwidth = 25)
}

#------------------------------------------------------------------------------
# Chimerism Calculation Functions
#------------------------------------------------------------------------------
#' Create Dataframe for Chimerism Calculation
#' @param GroupSelect Group identifier
#' @param AlleleTable Variant table (default: CommonAllelesGT)
#' @param DPcutoff Coverage threshold (default: 200)
#' @return Dataframe of informative sites
CreateChimerismDf <- function(GroupSelect, 
                            AlleleTable = CommonAllelesGT,
                            DPcutoff = 200) {
  # Find reference sample (first CD34neg donor)
  ComparisonIndiv <- AlleleTable %>% 
    filter(Group == GroupSelect & 
           SampleType == "DONOR" & 
           CellPopulation == "CD34neg") %>%
    arrange(Indiv) %>% 
    pull(Indiv) %>% 
    head(1)
  
  AlleleTable %>%
    filter(CHROM != "chrX") %>%
    filter(Group == GroupSelect) %>%
    group_by(Indiv) %>%
    filter(median(gt_DP) > DPcutoff) %>%
    group_by(ID) %>%
    mutate(ideal = case_when(
      gt_GT[Indiv == ComparisonIndiv] %in% c("0/0", "0|0") ~ 0,
      gt_GT[Indiv == ComparisonIndiv] %in% c("0/1", "0|1") ~ 0.5,
      gt_GT[Indiv == ComparisonIndiv] %in% c("1/1", "1|1") ~ 1)) %>%
    filter(ideal != 0.5) %>%
    mutate(diffIdeal = abs(ideal-VAF)) %>%
    filter(!grepl("Donor", SampleTimepoint)) %>%
    filter(!all(diffIdeal < 0.002))
}

#' Write VCF of Chimeric Sites
#' @param ChimerismDf Chimerism dataframe
#' @param path Output path
#' @param vcf VCF object (default: germlinevcf)
WriteChimericSitesVcf <- function(ChimerismDf, path, vcf = germlinevcf) {
  GroupSelect <- ChimerismDf$Group[1]
  sites <- ChimerismDf %>%
    pull(ID) %>% 
    unique()
  
  sitesvcf <- vcf[vcf@fix[,3] %in% sites]
  sitesvcf <- sitesvcf[samples = 1]
  write.vcf(sitesvcf, paste0(path, GroupSelect, ".ChimericSites.vcf.gz"))
}

#' Read Force Called Results
#' @param vcfpath Path to force called VCF
#' @param ChimerismDf Original chimerism dataframe
#' @return Processed force called results
ReadChimerismRecall <- function(vcfpath, ChimerismDf) {
  vcf <- read.vcfR(vcfpath) %>%
    vcfR2tidy(format_fields = c("GT", "AD", "DP"))
  
  recall <- vcf$fix[,1:6] %>%
    full_join(vcf$gt) %>% 
    filter(!grepl(",", ALT)) %>%
    mutate(VAF = VAFfromAD(gt_AD)) %>%
    left_join(ChimerismDf %>% select(ID, ideal) %>% unique()) %>%
    mutate(diffIdeal = abs(ideal-VAF)) %>%
    left_join(MetadataIDs, by = c("Indiv" = "Library_ID"))
  
  return(recall)
}

#' Cluster Chimeric Sites
#' @param recall Recalled variants
#' @param ChimerismDf Original chimerism data
#' @param TCellsOnly Only use T cell samples
#' @param k Number of clusters
#' @param DPcutoff Coverage threshold
#' @return Clustered sites
ClusterChimericSites <- function(recall, ChimerismDf, 
                               TCellsOnly = FALSE, k = 3,
                               DPcutoff = 100) {
  GroupSelect <- ChimerismDf$Group[1]
  
  samplestouse <- ChimerismDf %>% 
    ungroup %>% 
    select(Indiv, CellPopulation) %>% 
    unique() %>% 
    pull(Indiv, CellPopulation) 
  
  if(TCellsOnly) {
    samplestouse <- samplestouse[names(samplestouse) == "T"]
  }
  
  clustdf <- recall %>%
    filter(Indiv %in% samplestouse) %>%
    filter(gt_DP > DPcutoff) %>%
    select(ID, Indiv, diffIdeal) %>%
    pivot_wider(names_from = "Indiv", values_from = "diffIdeal") %>%
    na.omit()
  
  kclust <- clustdf %>%
    select(-ID) %>%
    kmeans(centers = k)
  
  clustdf <- augment(kclust, clustdf)
  
  p <- clustdf %>%
    ggpairs(columns = 2:(ncol(clustdf)-1),
            upper = "blank",
            aes(color = .cluster),
            progress = FALSE)
  
  ggsave(paste0("Figs/calculatechimerism/",
                GroupSelect,
                ".ggpairs.png"),
         plot = p,
         height = (ncol(clustdf)-2)*2,
         width = (ncol(clustdf)-2)*2)
  
  return(clustdf)
}

#' Calculate Final Chimerism Values
#' @param Clusters Clustered sites
#' @param recall Recalled variants
#' @param cluster Cluster to use
#' @param DPcutoff Coverage threshold
#' @return Chimerism calculations
CalculateChimerism <- function(Clusters, recall,
                             cluster, DPcutoff = 500) {
  GroupSelect <- recall$Group[1]
  IDs <- Clusters %>% 
    filter(.cluster == cluster) %>%
    pull(ID)
  
  ChimerismALL <- recall %>% 
    filter(ID %in% IDs) %>%
    na.omit() %>%
    group_by(Indiv) %>%
    summarise(ALL_Chimerism = mean(diffIdeal*2),
              ALL_SD = sd(diffIdeal*2),
              ALL_n = n(),
              ALL_meanDP = mean(gt_DP))
  
  recall %>% 
    filter(ID %in% IDs) %>%
    filter(gt_DP > DPcutoff) %>%
    group_by(Indiv) %>%
    summarise(Chimerism = mean(diffIdeal*2),
              SD = sd(diffIdeal*2),
              n = n(),
              meanDP = mean(gt_DP)) %>%
    full_join(ChimerismALL) %>%
    mutate(Group = GroupSelect)
}

#------------------------------------------------------------------------------
# Example Usage for One Group
#------------------------------------------------------------------------------
# Check for sample swaps in group
SampleSwapCheck("DR01")

# Create chimerism dataframe
chim_DR01 <- CreateChimerismDf("DR01")

# Write sites for force calling
WriteChimericSitesVcf(chim_DR01, "path/to/output/")

# NOTE: At this point, run FCchimerism.sh to force call sites
# Key command from shell script:
# gatk Mutect2 \
#   -R $genome_file \
#   -L DR01.ChimericSites.vcf.gz \
#   -alleles DR01.ChimericSites.vcf.gz \
#   -O DR01.ChimericSites.recall.vcf \
#   -I sample1.bam -I sample2.bam [...]

# gatk VariantAnnotator \
#   -V DR01.ChimericSites.recall.vcf \
#   -O DR01.ChimericSites.recall.ID.vcf \
#   -D $dbsnp_anno

# Read force called results
recall_DR01 <- ReadChimerismRecall("DR01.ChimericSites.recall.vcf", 
                                      chim_DR01)

# Cluster sites and visualize
clusters_DR01 <- ClusterChimericSites(recall_DR01, chim_DR01)

# Calculate final chimerism
# Choose cluster based on plots
chimerism_DR01 <- CalculateChimerism(clusters_DR01, 
                                     recall_DR01, 
                                     cluster = 2)

# View results
print(chimerism_DR01)