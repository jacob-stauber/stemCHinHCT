#!/bin/bash

##########################################################################################
##
## config.sh
##
##########################################################################################

trimmomatic_adapters=/gs/gsfs0/users/jstaube1/miniconda3/envs/TSA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
donor_pon=/gs/gsfs0/users/jstaube1/CHinHCT/RESULTS/FINAL/PON_20231128/pon_donor.vcf.gz
bashrc_path=/gs/gsfs0/users/jstaube1/.bashrc

genome_dir=/gs/gsfs0/users/jstaube1/CHinHCT/ref
genome_file=$genome_dir/hg38_align.fa
genomeindex=$genome_dir/hg38_align
dbsnp_file=$genome_dir/dbsnp_146.hg38.vcf
KGsnp_file=$genome_dir/1000G_phase1.snps.high_confidence.hg38.vcf
hapmap_file=$genome_dir/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
omni_file=$genome_dir/1000G_omni2.5.hg38.vcf.gz
indels_file=$genome_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf
axiom_file=$genome_dir/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
interval_file_bed=$genome_dir/targets_panelA.hg38.sort.bed
interval_file=$genome_dir/targets_panelA.hg38.sort.interval_list
refseq_file=$genome_dir/hg38.PanelATargetGenes.MANE.refseq
pon=$genome_dir/1000g_pon.hg38.vcf.gz
germline_resource=$genome_dir/af-only-gnomad.hg38.vcf.gz
common_snps=$genome_dir/small_exac_common_3.hg38.vcf.gz
funcotator_somatic=$genome_dir/funcotator_dataSources.v1.7.20200521s/
funcotator_germline=$genome_dir/funcotator_dataSources.v1.7.20200521g/
dbsnp_anno=$funcotator_somatic/dbsnp/hg38/hg38_All_20180418.vcf.gz

