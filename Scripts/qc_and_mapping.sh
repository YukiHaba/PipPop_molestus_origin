#!/bin/bash

##############################
## 0. Quality Control (FastQC)
##############################
# Run FastQC on all paired-end FASTQ files to assess initial read quality.
fastqc *.fq.gz -t 2 -o FastQCfiles 

##############################
## 1. Read Trimming (Trimmomatic)
##############################
# Define the adapter sequence file (e.g., Nextera adapters).
ADAPTER_SEQ_FILE="/path/to/Trimmomatic/adapters/NexteraPE-PE.fa"

# Run Trimmomatic in paired-end mode with 5 threads.
# Input: raw FASTQ files for read 1 and read 2.
# Output: paired and unpaired files for both forward and reverse reads.
java -jar /path/to/trimmomatic.jar \
    PE -threads 5 \
    ${SAMPLE_NAME}_R1.fastq.gz \           # Input FASTQ for read 1 (replace naming convention as needed)
    ${SAMPLE_NAME}_R2.fastq.gz \           # Input FASTQ for read 2
    ${SAMPLE_NAME}_forward_paired.fq.gz \  # Output: forward paired reads
    ${SAMPLE_NAME}_forward_unpaired.fq.gz \# Output: forward unpaired reads
    ${SAMPLE_NAME}_reverse_paired.fq.gz \  # Output: reverse paired reads
    ${SAMPLE_NAME}_reverse_unpaired.fq.gz \# Output: reverse unpaired reads
    ILLUMINACLIP:${ADAPTER_SEQ_FILE}:2:30:10 \  # Adapter clipping with specified parameters
    MINLEN:36                             # Discard reads shorter than 36 bases

##############################
## 2. Post-Trimming Quality Control (FastQC)
##############################
# Re-run FastQC on trimmed paired files to check quality improvements.
fastqc *.fq.gz -t 2 -o FastQCfiles

##############################
## 3. Read Alignment (BWA-MEM)
##############################
# Define folder containing input FASTQ files and/or output BAM files.
DATA_FOLDER="/path/to/data_folder"

# Specify the reference genome file.
REFERENCE="/path/to/reference_genome.fna"

sbatch -J "bwa.mem ${SAMPLE_NAME}" \
	/path/to/bwa_mem_batch_script.sh \
	${SAMPLE_NAME} \
	${DATA_FOLDER}/${SAMPLE_NAME}_forward_paired.fq.gz \    # Forward paired reads
	${DATA_FOLDER}/${SAMPLE_NAME}_reverse_paired.fq.gz \    # Reverse paired reads
	${DATA_FOLDER}/${SAMPLE_NAME}_forward_unpaired.fq.gz \  # Forward unpaired reads
	${DATA_FOLDER}/${SAMPLE_NAME}_reverse_unpaired.fq.gz \  # Reverse unpaired reads
	${REFERENCE}                                          # Reference genome

##############################
## 4. Marking Duplicates (Picard)
##############################
# Create directory to store duplication metrics.
mkdir -p MarkDupMetricsFile

# Run Picard MarkDuplicates on the sorted BAM file.
java -jar /path/to/picard.jar \
    MarkDuplicates \
    I=${SAMPLE_NAME}.sorted.bam \                      # Input sorted BAM file (assumes prior sorting)
    O=/tmp/${SAMPLE_NAME}.sorted.markdup.bam \         # Output BAM with duplicates marked
    M=/tmp/${SAMPLE_NAME}.sorted.markdup.metrics.txt   # Metrics file for duplicate marking

# Index the duplicate-marked BAM file for downstream analyses.
samtools index /tmp/${SAMPLE_NAME}.sorted.markdup.bam

# Optionally, generate mapping statistics using samtools flagstat.
mkdir -p FlagstatFile
samtools flagstat /tmp/${SAMPLE_NAME}.sorted.markdup.bam > /tmp/${SAMPLE_NAME}.sorted.markdup.bam.flagstat

##############################
## 5. Indel Realignment (GATK)
##############################
# Define the input BAM file variable for this sample.
SAMPLE_BAM=${SAMPLE_NAME}.sorted  # Adjust naming if needed

# Step 5a: Identify target intervals for indel realignment.
java -Xmx20g -jar /path/to/GenomeAnalysisTK.jar \
    -nt 1 \
    -T RealignerTargetCreator \
    -R ${REFERENCE} \
    -I ${SAMPLE_BAM}.bam \                 # Input BAM file for the sample
    -o ${SAMPLE_BAM}.indel.list            # Output list of indel targets

# Step 5b: Perform local realignment around the identified indels.
java -Xmx20g -jar /path/to/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${REFERENCE} \
    -I ${SAMPLE_BAM}.bam \                 # Same input BAM file
    -targetIntervals ${SAMPLE_BAM}.indel.list \  # Use the indel intervals from previous step
    -o ${SAMPLE_BAM}.realigned.bam         # Output realigned BAM file

##############################
## 6. Coverage Calculation (mosdepth)
##############################
# Compute coverage statistics from the final realigned BAM file.
mosdepth \
    /tmp/${SAMPLE_NAME} \                    # Prefix for mosdepth output files
    ${SAMPLE_BAM}.realigned.bam \            # Input BAM file (post realignment)
    --d4                                     # Option for per-base depth (adjust options as needed)