#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=VQSR
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=8				# Run on N CPUs
#SBATCH --mem=36gb						# Job memory request
#SBATCH --time=24:00:00               	# Time limit hrs:min:sec
#SBATCH --output=%x_%j.log   	# Standard output and error log
pwd; hostname; date
echo "cores" $SLURM_CPUS_PER_TASK
echo "job" $SLURM_JOB_NAME
echo "ID" $SLURM_JOB_ID

HOME=/gs/gsfs0/users/jstaube1

echo '---- Preping Enviroment ----'
echo -e "$(date) \t timestamp: $(date +%s)" 

source $HOME/.bashrc
echo $PATH
echo $ID
echo $Fastq1

source ~/CHinHCT/scripts/config.hpc3.sh

conda activate gatk

echo '---- Joint Genotyping ----'
echo -e "$(date) \t timestamp: $(date +%s)" 

gatk --java-options "-Xmx24g -Xms4g" GenomicsDBImport \
	--genomicsdb-workspace-path CiH_ALL_db \
	--batch-size 50 \
	-L $interval_file \
	-ip 100 \
	--merge-input-intervals true \
	--sample-name-map CiH_ALL.sample_map \
	--tmp-dir tmp \
	--reader-threads 8
	
gatk --java-options "-Xmx24g" GenotypeGVCFs \
	-R $genome_file \
	-V gendb://CiH_ALL_db \
	-O CiH_ALL.vcf.gz
	
gatk MakeSitesOnlyVcf \
	-I CiH_ALL.vcf.gz \
	-O CiH_ALL_sitesonly.vcf.gz
	
echo '---- Calculating Variant Recalibration ----'
echo -e "$(date) \t timestamp: $(date +%s)" 
	
gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
    -V CiH_ALL_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR \
    -mode INDEL \
    --max-gaussians 1 \
    -resource:mills,known=false,training=true,truth=true,prior=12 $indels_file \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 $axiom_file \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $dbsnp_file \
    -O CiH_ALL.indels.recal \
    --tranches-file CiH_ALL.indels.tranches \
    --rscript-file CiH_ALL.indels.plots.R
    
gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
    -V CiH_ALL_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    --max-gaussians 1 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_file \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni_file \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $KGsnp_file \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_file \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
    -O CiH_ALL.snps.recal \
    --tranches-file CiH_ALL.snps.tranches \
    --rscript-file CiH_ALL.snps.plots.R

echo '---- Applying Recalibration ----' 
echo -e "$(date) \t timestamp: $(date +%s)" 
    
gatk --java-options "-Xmx24g -Xms5g" ApplyVQSR \
    -V CiH_ALL.vcf.gz \
    --recal-file CiH_ALL.indels.recal \
    --tranches-file CiH_ALL.indels.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL \
    -O CiH_ALL.indel.recalibrated.vcf.gz
    
gatk --java-options "-Xmx24g -Xms5g" ApplyVQSR \
    -V CiH_ALL.indel.recalibrated.vcf.gz \
    --recal-file CiH_ALL.snps.recal \
    --tranches-file CiH_ALL.snps.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode SNP \
    -O CiH_ALL.recalibrated.vcf.gz
    
echo '---- Variant Annotation ----'
echo -e "$(date) \t timestamp: $(date +%s)" 
    
gatk VariantAnnotator \
	-V CiH_ALL.recalibrated.vcf.gz \
	-O CiH_ALL.recalibrated.dbsnp.vcf.gz \
	-D $dbsnp_anno
	
gatk --java-options "-Xmx24G" Funcotator \
   -R $genome_file \
   -V CiH_ALL.recalibrated.dbsnp.vcf.gz  \
   -O CiH_ALL.recalibrated.funcotator.maf \
   --output-file-format MAF \
   --data-sources-path $funcotator_germline \
   --ref-version hg38
   
conda activate vcf2maf

gunzip -kd CiH_ALL.recalibrated.dbsnp.vcf.gz
	
vcf2maf.pl \
	--input-vcf CiH_ALL.recalibrated.dbsnp.vcf  \
	--output-maf CiH_ALL.recalibrated.vcf2maf.maf \
	--ref-fasta $genome_file \
	--cache-version 105

echo '---- Analysis Complete ----'
echo -e "$(date) \t timestamp: $(date +%s)" 


#https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
#https://gatk.broadinstitute.org/hc/en-us/articles/360035890831
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531612?id=1259
#https://gatk.broadinstitute.org/hc/en-us/articles/4402736812443-Which-training-sets-arguments-should-I-use-for-running-VQSR-

	


	
