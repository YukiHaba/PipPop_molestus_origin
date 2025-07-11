# Population Genomics Analysis Pipeline Scripts

This directory contains a sequential pipeline of bash scripts for population genomics analysis of NGS data. The scripts are designed for analysis of *Culex* mosquito genomic data but can be adapted for other organisms.

## Pipeline Overview

The pipeline consists of 5 scripts that should be run in numerical order with specific file dependencies:

### 0. Quality Control and Read Mapping (`0. qc_and_mapping.sh`)
- **Purpose**: Process raw sequencing reads through quality control and alignment
- **Steps**: FastQC → Trimmomatic → BWA-MEM → Picard MarkDuplicates → GATK IndelRealigner → mosdepth
- **Input**: Raw FASTQ files
- **Output**: `${SAMPLE_NAME}.realigned.bam` files with quality metrics

### 1. Variant Calling and Filtering (`1. variant_calling.sh`)
- **Purpose**: Call variants from aligned reads and apply quality filters
- **Steps**: bcftools mpileup/call → Quality filtering → Repeat masking → Accessible site extraction → LD pruning
- **Input**: Realigned BAM files from script 0
- **Output**: 
  - `good_biallelic_snps.rm.combined.accessible.vcf.gz` (main variant set)
  - `good_biallelic_snps.rm.combined.accessible.LDpruned.vcf.gz` (LD-pruned variants)

### 2. Kinship Analysis (`2. kinship_ngsRelate.sh`)
- **Purpose**: Calculate kinship coefficients between individuals
- **Steps**: ngsRelate analysis
- **Input**: `good_biallelic_snps.rm.combined.accessible.vcf.gz` from script 1
- **Output**: Relatedness estimates

### 3. Haplotype Phasing (`3. phasing.sh`)
- **Purpose**: Phase variants into haplotypes using read-based and statistical methods
- **Steps**: extractHAIRS → HAPCUT2 → whatshap statistics → SHAPEIT4
- **Input**: 
  - `good_biallelic_snps.rm.combined.accessible.vcf.gz` from script 1
  - `${SAMPLE_NAME}.realigned.bam` from script 0
- **Output**: `${SAMPLE_ID}.${CHR}.phased.bcf` (phased haplotypes)

### 4. Population Genomics Analyses (`4. popgen_analyses.sh`)
- **Purpose**: Population genetic analyses and demographic inference
- **Steps**: PCA → f3/f4 statistics → Nucleotide diversity → MSMC2 demographic history
- **Input**: 
  - `good_biallelic_snps.rm.combined.accessible.LDpruned.vcf.gz` from script 1
  - `${SAMPLE_ID}.${CHR}.phased.bcf` from script 3 (for MSMC2)
  - `${SAMPLE_NAME}.realigned.bam` from script 0 (for MSMC2)
- **Output**: Population structure and demographic history results

## Data

### Raw Sequencing Data
- Paired-end Illumina sequencing data (FASTQ files) - available from NCBI project PRJNA1209100

### Reference Files (provided in this repository)
- Reference genome: `GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna`
- Repeat regions: `GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_rm.out.gz`
- GTF annotation: `GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.gtf.gz`

### Additional Required Files
- Accessible genomic regions (BED format)
- Population assignment files
- Genetic map files (for phasing)

## Notes on Configuration

Each script contains a configuration section at the top with variables that should be modified for your system:

- File paths (input/output directories, reference genome, etc.)
- Analysis parameters (threads, quality thresholds, etc.)
- Tool paths and conda environments

## Basic Usage

1. Download raw FASTQ files from NCBI project PRJNA1209100
2. Modify configuration variables in each script
3. Use scripts sequentially: `0 → 1 → 2 → 3 → 4`
4. Monitor intermediate outputs and adjust parameters as needed

## Author

Yuki Haba
