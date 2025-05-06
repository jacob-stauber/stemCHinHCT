# Dynamics of Donor Stem Cell Clonal Hematopoiesis following Hematopoietic Cell Transplantation

**Authors:** Jacob Stauber, Oliver Bohorquez, Catherine Jankovic, Caner Saygin, Sriram Sundaravel, Jong Jeong, Satyajit Kosuri, Michael Bishop, Amittha Wickrema, John Greally, Ulrich Steidl

## Study Overview

This repository contains the data and code supporting the manuscript "Dynamics of Donor Stem Cell Clonal Hematopoiesis following Hematopoietic Cell Transplantation."

The study investigates how clonal hematopoiesis (CH) variants in allogeneic hematopoietic cell transplant (HCT) donors behave after transplantation. We used targeted deep sequencing to track CH clones in mobilized donor peripheral blood (fractionated into CD34+ HSPC-enriched and CD34- populations) and in matched longitudinal recipient samples (sorted into specific lineages and HSPCs) during the first year post-transplant.

## Code Description

The code provided in this repository was used for the bioinformatics analysis of targeted sequencing data. The overall workflow includes the following main stages:

1.  **Raw Sequencing Data Processing (`scripts/shell/TargetedSeqAnalysis.v3.sh`):**
    * Adapter trimming of raw FASTQ files.
    * Alignment to the GRCh38 reference genome.
    * UMI processing to generate consensus reads for error correction.
    * Initial somatic variant calling on a per-sample basis using GATK Mutect2.

2.  **Germline Variant Calling and Chimerism Calculation (`scripts/shell/JointCalling.sh`, `scripts/R/germline-chimerism.r`):**
    * Joint calling of germline variants across all samples using GATK HaplotypeCaller and GenomicsDBImport.
    * Variant Quality Score Recalibration (VQSR) to filter high-confidence germline SNPs.
    * Calculation of donor chimerism in recipient samples using informative germline SNPs.

3.  **Somatic Variant Refinement and Analysis (`scripts/R/01.prelim-analysis.r`, `scripts/R/02.forcecall.r`, `scripts/shell/ForceCall.sh`, `scripts/R/03.final-vars.r`):**
    * **Preliminary Filtering (`01.prelim-analysis.r`):** Initial filtering of somatic variant calls from `TargetedSeqAnalysis.v3.sh` based on read depth, allele frequency, and quality metrics.
    * **Force Calling Preparation (`02.forcecall.r`):** Identification of variants for deeper investigation across all samples within specific donor-recipient groups.
    * **Force Calling Execution (`ForceCall.sh`):** Re-running GATK Mutect2 in "force-call" mode to assess the presence and VAF of selected variants in every sample of a group, even if not initially called.
    * **Final Variant Curation (`03.final-vars.r`):** Integration of force-called results, chimerism data, and annotations. Application of final filters and manual review (based on an external review file) to produce the definitive list of CH variants for the study.

The `scripts/shell/config.hpc3.sh` file contains paths and settings used by the shell scripts, and `environment/TSA_condaenv.yml` defines the Conda environment for the primary analysis pipeline.

## Data File: `FinalVars_nonsynonymous.csv`

This file contains the final dataset of all donor derived non-synonymous somatic variants tracked in the study, after all filtering and processing steps. The data is provided in 'tidy' format.


### Key Columns:

* `Group`: Donor-recipient pair ID.
* `SampleType`: Indicates if the sample is from the "DONOR" or recipient "BLOOD" or "MARROW".
* `SampleTimepoint`: Descriptive timepoint of sample collection.
* `Time`: Numeric time relative to transplant in days.
* `CellPopulation`: The cell population from which DNA was extracted.
* `Library_ID`: Unique identifier for the sequencing library.
* `VCFVAR`: Variant representation in CHR-POS-REF-ALT format.
* `Hugo_Symbol`
* `HGVSc`: Coding sequence change.
* `HGVSp_Short`: Protein sequence change.
* `t_depth`: Total read depth in the specified sample.
* `t_alt_count`: Number of alt reads in the specified sample.
* `VAF`: Variant Allele Frequency.
* `Variant_Classification`: Functional classification (e.g., "Missense_Mutation", "Frame_Shift_Del", "Nonsense_Mutation", "Splice_Site").
* `Variant_Type`: Type of sequence alteration (e.g., "SNP", "DEL", "INS").
* `Existing_variation`: Known identifiers for the variant from public databases (e.g., dbSNP rsID, COSMIC ID).
* `OrigFilter`: The filter status from the initial, per-sample Mutect2 call.
* `IN_SCREEN`: Boolean flag indicating if the variant was identified during the screening phase.
* `DonorChimerism`: The calculated donor chimerism for the recipient sample at the given timepoint and cell population.
* `ChimerismSD`: The standard deviation of the chimerism measurement.