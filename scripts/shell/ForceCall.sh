#!/bin/bash
#SBATCH -p normal                        # Batch partition: test, quick, ht, normal, unlimited
#SBATCH --job-name=jupyter              # Job name
#SBATCH --ntasks=1                      # Number of task to run in parallel
#SBATCH --cpus-per-task=16               # Run on N CPUs for multithreading
#SBATCH --mem=72gb                      # Job memory request
#SBATCH --time=48:00:00                 # Time limit hrs:min:sec
#SBATCH --output=cd3FC_%j.log         # Standard output and error log - %j will add job ID to file name
pwd; hostname; date

source ~/CHinHCT/scripts/config.hpc3.sh
source $bashrc_path
RAM=64


while read -r LINE; do
	while IFS=',' read -ra VARS; do
		GROUP="${VARS[0]}"
		INPUT="${VARS[1]}"
		NAME="$(basename $INPUT .bam)"
		echo $GROUP
		echo "$INPUT"
		echo $NAME
		
		#convert site MAF exported from R into vcf appropriate for Mutect2 input
		#skip if already done per group
		if [ ! -f $GROUP/$GROUP.sites.sort.vcf.gz ]; then
			mkdir $GROUP
		
			conda activate vcf2maf
			maf2vcf.pl \
				--input-maf ./inMAFs/$GROUP.sites.maf \
				--output-dir $GROUP --ref-fasta $genome_file
		
			conda activate TSA
		
			vcf-sort $GROUP/$GROUP.sites.vcf > $GROUP/$GROUP.sites.sort.vcf
		
			bgzip $GROUP/$GROUP.sites.sort.vcf
			tabix $GROUP/$GROUP.sites.sort.vcf.gz
		fi
			
		#force call by sample and annotate
		#skip filter
		conda activate TSA
		
		gatk --java-options "-Xmx${RAM}G" Mutect2 \
			-R $genome_file \
			-L $GROUP/$GROUP.sites.sort.vcf.gz \
			-alleles $GROUP/$GROUP.sites.sort.vcf.gz \
			-O $GROUP/$GROUP.$NAME.FC.vcf \
			-I $INPUT \
			--germline-resource $germline_resource
			
		gatk LeftAlignAndTrimVariants \
			-R $genome_file \
			-V $GROUP/$GROUP.$NAME.FC.vcf \
			-O $GROUP/$GROUP.$NAME.FC.split.vcf \
			--split-multi-allelics
		
		conda activate vcf2maf
		vcf2maf.pl \
			--input-vcf $GROUP/$GROUP.$NAME.FC.split.vcf \
			--output-maf $GROUP/$GROUP.$NAME.FC.split.MAF \
			--ref-fasta $genome_file \
			--cache-version 105 \
			--tumor-id $NAME	
		
		MERGE_FILE=$GROUP/$GROUP.FC.cD3.MAF
		#merge all MAFs per group	
		if [ ! -f $MERGE_FILE ]; then
			head -2 $GROUP/$GROUP.$NAME.FC.split.MAF > $MERGE_FILE
		fi
		tail -n +3 -q $GROUP/$GROUP.$NAME.FC.split.MAF >> $MERGE_FILE

		echo
		echo
	done <<< "$LINE"
done < ./filelist_cD3.indiv.missing.txt
