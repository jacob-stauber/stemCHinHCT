#!/bin/bash

##########################################################################################
##
## TargetedSeqAnalysis.sh
##
## Main workflow
## For use on HPC3.0 - Jacob Stauber
## Last edited 02/01/2024
## Adds Support for multiple fastq and calculates start/stop dup rate
##
##########################################################################################

usage()
{
	echo "  Usage: $0 "
	echo "	-n, --name               Name of the sample."
	echo "	-D,    					 Path to fastq directory. Multiple directories can be provided seperated by ';'. Files must end in _1.fq.gz,_2.fq.gz,_3.fq.gz"
	echo "	-cD, --consensusdepth    Consensus Depth passed to FilterConsensusReads"
	echo "	-c, --config             Path to configuration file. Optional."
	echo "	-t, --threads            Number of CPU threads. Optional. Defaults to 8."
	echo "	-r, --RAM                Amount of Gb RAM. Optional. Defaults to 32."
	echo "	-temp, --temp_dir        Path to temporary directory."
	echo "	-filt, --filtering       Set to 'all' (AF >= 0.05, , Variant in Tumor >= 2, Variant in Normal <= 1, Coverage >= 5), 'hard' (AF >= 0.1, Variant in Tumor >= 3, Variant in Normal = 0, Coverage >= 10) or 'none' (no filters). Optional. Defaults to 'hard'."
	echo "	--help                   Show this help."
  exit 1
}

# default parameters
fastq_1=
fastq_2=
fastq_umi=
bam_normal=
threads=8
RAM=32
temp_dir=
filtering=hard

# parse parameters
if [ "$1" = "" ]; then usage; fi
while [ "$1" != "" ]; do case $1 in
	-n|--name) shift;name="$1";;
	-D) shift;DIR="$1";;
	-cD|--consensusdepth) shift;ConsensusDepth="$1";;
	-c|--config) shift;config_file="$1";;
	-t|--threads) shift;threads="$1";;
	-r|--RAM) shift;RAM="$1";;
	-temp|--temp_dir) shift;temp_dir="$1";;
	-filt|--filtering) shift;filtering="$1";;
    --help) usage;shift;;
	*) usage;shift;;
esac; shift; done


#reading configuration from $config_file
source $config_file
source $bashrc_path #fixes potential issues with conda

echo '---- Starting Targeted Seq Analysis ----'
echo -e "$(date) \t timestamp: $(date +%s)"

echo '---- Creating directories ----'
echo -e "$(date) \t timestamp: $(date +%s)"
mkdir -p $name/
mkdir -p $name/pipeline
mkdir -p $name/fastq
mkdir -p $name/results/QC
mkdir -p $name/results/bam
mkdir -p $name/results/Mutect2
mkdir -p $name/results/HaplotypeCaller

echo '---- Parsing Inputs ----'
echo -e "$(date) \t timestamp: $(date +%s)"
IFS=';' read -ra ARR <<< "$DIR"
fastq_1=$(for i in ${ARR[@]}; do echo $i/*_1.fq.gz; done)
fastq_2=$(for i in ${ARR[@]}; do echo $i/*_3.fq.gz; done)
fastq_umi=$(for i in ${ARR[@]}; do echo $i/*_2.fq.gz; done)


MAX_RECORDS_IN_RAM=$(expr $RAM \* 250000)
HASH_TABLE_SIZE=$((RAM*1000000000/500))

echo '---- Starting UMI Targeted Seq Analysis ----' | tee -a $name/results/QC/$name.report.txt
echo Starting pipeline using these settings: | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
echo Running sample named $name | tee -a $name/results/QC/$name.report.txt
echo Using $fastq_1 and $fastq_2  | tee -a $name/results/QC/$name.report.txt
echo Using $fastq_umi for UMI fastqs | tee -a $name/results/QC/$name.report.txt
echo Reading configuration file from $config_file | tee -a $name/results/QC/$name.report.txt
#echo Setting location of repository to $repository_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location of genome to $genome_dir | tee -a $name/results/QC/$name.report.txt
echo Setting location for temporary files to $temp_dir| tee -a $name/results/QC/$name.report.txt
echo Setting cD threshold to $ConsensusDepth| tee -a $name/results/QC/$name.report.txt
echo $filtering is setting for filtering of SNV calls | tee -a $name/results/QC/$name.report.txt
echo
echo Starting workflow using $threads CPU-threads and $RAM GB of RAM | tee -a $name/results/QC/$name.report.txt

#rerouting STDERR to report file
exec 2>> $name/results/QC/$name.report.txt



echo '---- Creating directories ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

echo '---- Loading Conda Enviroment ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

conda activate TSA

cp $config_file $name/pipeline/
cp "$0" $name/pipeline/ #copy script
 
echo '---- Copying raw data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

cat $fastq_1 > $name/fastq/$name.R1.fastq.gz
cat $fastq_2 > $name/fastq/$name.R2.fastq.gz
cat $fastq_umi > $name/fastq/$name.UMI.fastq.gz

echo '---- Running FastQC before trimming ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

fastqc -t $threads \
	$name/fastq/$name.R1.fastq.gz \
	$name/fastq/$name.R2.fastq.gz \
	--outdir=$name/results/QC

echo '---- Trimming reads ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

trimmomatic PE \
	-threads $threads -phred33 \
	$name/fastq/$name.R1.fastq.gz \
	$name/fastq/$name.R2.fastq.gz \
	$name/fastq/$name.R1.passed.fastq.gz \
	$name/fastq/$name.R1.not_passed.fastq.gz \
	$name/fastq/$name.R2.passed.fastq.gz \
	$name/fastq/$name.R2.not_passed.fastq.gz \
	ILLUMINACLIP:$trimmomatic_adapters:2:30:10:8:TRUE

echo '---- Running FastQC after trimming ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

fastqc -t $threads \
	$name/fastq/$name.R1.passed.fastq.gz \
	$name/fastq/$name.R2.passed.fastq.gz \
	--outdir=$name/results/QC
		
echo '---- Mapping trimmed reads ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

bwa mem -t $threads $genomeindex -v 1 \
	$name/fastq/$name.R1.passed.fastq.gz \
	$name/fastq/$name.R2.passed.fastq.gz \
	| samtools view -o $name/results/bam/$name.mapped.bam -
	
gatk --java-options "-Xmx${RAM}G" FixMateInformation \
	-I $name/results/bam/$name.mapped.bam  \
	-O $name/results/bam/$name.mapped.fixed.bam

gatk --java-options "-Xmx${RAM}G" AddOrReplaceReadGroups \
	-I $name/results/bam/$name.mapped.fixed.bam  \
	-O $name/results/bam/$name.mapped.fixed.addRG.bam \
	-LB $name -PL ILLUMINA -PU $name -SM $name \
	--SORT_ORDER coordinate --CREATE_INDEX true
		
	
echo '---- Removing fastq files ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

rm $name/fastq/$name.R1.passed.fastq.gz
rm $name/fastq/$name.R1.not_passed.fastq.gz
rm $name/fastq/$name.R2.passed.fastq.gz
rm $name/fastq/$name.R2.not_passed.fastq.gz

echo '---- Base recalibration ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

gatk --java-options "-Xmx${RAM}G" BaseRecalibrator \
	-R $genome_file \
	-I $name/results/bam/$name.mapped.fixed.addRG.bam	 \
	--known-sites $dbsnp_file \
	--known-sites $KGsnp_file \
	--known-sites $indels_file \
	-L $interval_file \
	-ip 100 \
	-O $name/results/QC/$name.pre.recal.table

gatk --java-options "-Xmx${RAM}G" ApplyBQSR \
	-R $genome_file \
	-I $name/results/bam/$name.mapped.fixed.addRG.bam\
	-O $name/results/bam/$name.mapped.fixed.addRG.recal.bam \
	-bqsr $name/results/QC/$name.pre.recal.table

gatk --java-options "-Xmx${RAM}G" BaseRecalibrator \
	-R $genome_file \
	-I $name/results/bam/$name.mapped.fixed.addRG.recal.bam  \
	--known-sites $dbsnp_file \
	--known-sites $KGsnp_file \
	--known-sites $indels_file \
	-L $interval_file \
	-ip 100 \
	-O $name/results/QC/$name.post.recal.table	
	
gatk --java-options "-Xmx${RAM}G" AnalyzeCovariates \
    -before $name/results/QC/$name.pre.recal.table \
    -after $name/results/QC/$name.post.recal.table \
    -plots $name/results/QC/$name.AnalyzeCovariates.pdf \
    -csv $name/results/QC/$name.AnalyzeCovariates.csv

echo '---- Adding UMIs to mapped reads ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

fgbio -Xmx${RAM}G -XX:+AggressiveHeap AnnotateBamWithUmis \
	-i $name/results/bam/$name.mapped.fixed.addRG.recal.bam  \
	-f $name/fastq/$name.UMI.fastq.gz \
	-o $name/results/bam/$name.mapped.fixed.addRG.recal.addUMI.bam 

echo '---- Grouping reads by UMI ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

fgbio -Xmx${RAM}G -XX:+AggressiveHeap GroupReadsByUmi \
	--input=$name/results/bam/$name.mapped.fixed.addRG.recal.addUMI.bam  \
	--output=$name/results/bam/$name.mapped.fixed.addRG.recal.addUMI.grouped.bam  \
	--family-size-histogram=$name/results/QC/$name.UMImetrics.txt \
	--strategy=edit \
	--edits=1 \
	--min-map-q=20

echo '---- Generating consesus reads ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

fgbio -Xmx${RAM}G -XX:+AggressiveHeap CallMolecularConsensusReads \
	--input=$name/results/bam/$name.mapped.fixed.addRG.recal.addUMI.grouped.bam \
	--output=$name/results/bam/$name.consensus.unmapped.bam \
	--sort-order=queryname \
	--error-rate-post-umi=30 \
	--min-reads=1

# java -Xmx${RAM}G -XX:+AggressiveHeap -jar $fgbio_jar FilterConsensusReads \
# 	--input=$name/results/bam/$name.consensus.unmapped.bam  \
# 	--output=$name/results/bam/$name.consensus.unmapped.filtered.bam \
# 	--ref=$genome_file \
# 	--reverse-per-base-tags=true \
# 	--min-reads=1 \
# 	--max-read-error-rate=0.05 \
# 	--min-base-quality=40 \
# 	--max-base-error-rate=0.1 \
# 	--max-no-call-fraction=0.2
	
	
fgbio -Xmx${RAM}G -XX:+AggressiveHeap FilterConsensusReads \
	--input=$name/results/bam/$name.consensus.unmapped.bam  \
	--output=$name/results/bam/$name.consensus.unmapped.filtered.bam \
	--ref=$genome_file \
	--reverse-per-base-tags=true \
	--min-reads=$ConsensusDepth \
	--max-read-error-rate=0.05 \
	--min-base-quality=20 \
	--max-base-error-rate=0.1 \
	--max-no-call-fraction=0.2
	
gatk --java-options "-Xmx${RAM}G" SortSam \
	-I $name/results/bam/$name.consensus.unmapped.filtered.bam \
	-O $name/results/bam/$name.consensus.unmapped.filtered.sorted.bam \
	--SORT_ORDER queryname	
	
echo '---- Remapping consensus reads ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

gatk --java-options "-Xmx${RAM}G" SamToFastq \
	-I $name/results/bam/$name.consensus.unmapped.filtered.sorted.bam \
	-F /dev/stdout \
	-INTER true \
		| bwa mem -p -t $threads -v 1 $genomeindex /dev/stdin \
		| gatk --java-options "-Xmx${RAM}G" MergeBamAlignment \
			-UNMAPPED $name/results/bam/$name.consensus.unmapped.filtered.sorted.bam \
			-ALIGNED /dev/stdin \
			-O $name/results/bam/$name.bam \
			-R $genome_file \
			-SO coordinate \
			--ALIGNER_PROPER_PAIR_FLAGS true \
			-MAX_GAPS -1 \
			-ORIENTATIONS FR \
			--VALIDATION_STRINGENCY SILENT \
			--CREATE_INDEX true


echo '---- Collecting Alignment Metrics ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

gatk --java-options "-Xmx${RAM}G" CollectMultipleMetrics \
	-R $genome_file \
	-I $name/results/bam/$name.bam \
	-O $name/results/QC/$name.bam.metrics \
	--PROGRAM CollectAlignmentSummaryMetrics \
	--PROGRAM CollectInsertSizeMetrics \
	--PROGRAM QualityScoreDistribution \
	--PROGRAM MeanQualityByCycle \
	--PROGRAM CollectBaseDistributionByCycle \
	--PROGRAM CollectGcBiasMetrics \
	--PROGRAM CollectSequencingArtifactMetrics \
	--PROGRAM CollectQualityYieldMetrics

gatk --java-options "-Xmx${RAM}G" CollectHsMetrics \
	-R $genome_file \
	-I $name/results/bam/$name.bam \
	-O $name/results/QC/$name.bam.HsMetrics \
	--BAIT_INTERVALS $interval_file \
	--TARGET_INTERVALS $interval_file

gatk --java-options "-Xmx${RAM}G" DepthOfCoverage \
	-R $genome_file \
	-I $name/results/bam/$name.bam \
	-O $name/results/QC/$name.bam.DoC \
	-L $interval_file \
	-gene-list $refseq_file

mosdepth -n --by $interval_file_bed \
	$name/results/QC/$name.bam \
	$name/results/bam/$name.bam


samtools idxstats $name/results/bam/$name.bam \
	> $name/results/QC/$name.bam.idxstats
	
samtools flagstat $name/results/bam/$name.bam \
	> $name/results/QC/$name.bam.flagstat

echo '---- Summarizing quality control data ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

multiqc $name/results/QC -n $name -o $name/results/QC/


echo '---- Running Mutect2 (single-sample) ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt


gatk --java-options "-Xmx${RAM}G" Mutect2 \
	--native-pair-hmm-threads $threads \
	-R $genome_file \
	-I $name/results/bam/$name.bam \
	-O $name/results/Mutect2/$name.m2.vcf \
	-L $interval_file \
	-ip 100 \
	-pon $donor_pon \
	--germline-resource $germline_resource \
	--f1r2-tar-gz $name/results/Mutect2/$name.f1r2.tar.gz \
	-bamout $name/results/Mutect2/$name.m2.bam
	
echo '---- Filtering Mutect2 Calls ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	
gatk --java-options "-Xmx${RAM}G" GetPileupSummaries \
	-I $name/results/bam/$name.bam \
	-V $common_snps \
	-L $interval_file \
	-ip 100 \
	-O $name/results/Mutect2/$name.pileups.table
	
gatk --java-options "-Xmx${RAM}G" CalculateContamination \
	-I $name/results/Mutect2/$name.pileups.table \
	-segments $name/results/Mutect2/$name.segments.table \
	-O $name/results/Mutect2/$name.contamination.table
	
gatk --java-options "-Xmx${RAM}G" LearnReadOrientationModel \
	-I $name/results/Mutect2/$name.f1r2.tar.gz\
	-O $name/results/Mutect2/$name.artifact-prior.tar.gz
	

gatk --java-options "-Xmx${RAM}G" FilterMutectCalls \
	-R $genome_file \
	-V $name/results/Mutect2/$name.m2.vcf \
	--max-n-ratio 0.05 \
	--contamination-table $name/results/Mutect2/$name.contamination.table \
	--orientation-bias-artifact-priors $name/results/Mutect2/$name.artifact-prior.tar.gz \
	-O $name/results/Mutect2/$name.m2.filt.vcf
	
vcftools \
	--vcf $name/results/Mutect2/$name.m2.filt.vcf \
	--remove-filtered-all --recode --stdout \
	> $name/results/Mutect2/$name.m2.filt.PASS.vcf
	
echo '---- Annotating Mutect2 Calls ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt
	
conda activate vcf2maf
PERL5LIB=""	
	
vcf2maf.pl \
	--input-vcf $name/results/Mutect2/$name.m2.filt.PASS.vcf \
	--output-maf $name/results/Mutect2/$name.m2.filt.PASS.maf \
	--tumor-id $name \
	--ref-fasta $genome_file \
	--cache-version 105 \
	--verbose

vcf2maf.pl \
	--input-vcf $name/results/Mutect2/$name.m2.filt.vcf \
	--output-maf $name/results/Mutect2/$name.m2.filt.maf \
	--tumor-id $name \
	--ref-fasta $genome_file \
	--cache-version 105 \
	--verbose
	
echo '---- Single-sample GVCF calling ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt	

conda activate gatk
   
gatk --java-options "-Xmx${RAM}G" HaplotypeCaller \
	--native-pair-hmm-threads $threads \
	-R $genome_file \
	-L $interval_file \
	-ip 100 \
	-I $name/results/bam/$name.bam \
    -O $name/results/HaplotypeCaller/$name.g.vcf.gz \
    -ERC GVCF

echo '---- Calculating Start/Stop Dup rate ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt	

conda activate samtools

samtools sort -n $name/results/bam/$name.bam > $name/results/bam/$name.ST.sort.bam
samtools fixmate -m $name/results/bam/$name.ST.sort.bam $name/results/bam/$name.ST.fix.bam
samtools sort $name/results/bam/$name.ST.fix.bam > $name/results/bam/$name.ST.resort.bam
samtools markdup --duplicate-count $name/results/bam/$name.ST.resort.bam $name/results/bam/$name.mkdup.bam
TAG=dc
samtools view $name/results/bam/$name.mkdup.bam | grep "$TAG:." | perl -pe 's/(^.+?)\t.*'$TAG':.:(.+?)/$1\t$2/g' > $name/results/QC/$name.SSduprate.tsv

awk '{print $2}' $name/results/QC/$name.SSduprate.tsv | sort -n | uniq -c > $name/results/QC/$name.SSduprate.hist.txt


#### For Panel of Normals ####
# java -Xmx${RAM}G -jar $GATK_jar Mutect2 \
# 	--native-pair-hmm-threads $threads \
# 	-R $genome_file \
# 	-I $name/results/bam/$name.bam \
# 	-O $name/results/Mutect2/$name.m2.vcf \
# 	-L $interval_file \
# 	-ip 100 \
# 	-max-mnp-distance 0
#### END PON Block ####


echo '---- Finished analysis of sample '$name' ----' | tee -a $name/results/QC/$name.report.txt
echo -e "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

wc -l $name/results/Mutect2/$name.m2.filt.vcf >> CallCount.txt #useful for checking progress

exit 0