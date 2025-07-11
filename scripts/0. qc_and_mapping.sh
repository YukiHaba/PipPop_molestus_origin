#!/bin/bash

################################################################################
# SCRIPT: 0. Quality Control and Read Mapping Pipeline
# AUTHOR: Yuki Haba
# REFERENCE: Yuki Haba et al. 2025 Ancient origin of an urban underground 
#            mosquito (PMID: 39975080)
# DESCRIPTION: Pipeline for NGS data quality control, trimming,
#              alignment, duplicate marking, indel realignment, and coverage
#              calculation. Designed for paired-end Illumina sequencing data.
################################################################################

# Exit on any error
set -e

#===============================================================================
# CONFIGURATION VARIABLES - MODIFY THESE PATHS FOR YOUR SYSTEM
#===============================================================================

# Sample information (required variables)
SAMPLE_NAME="${SAMPLE_NAME:-sample_name}"  # Sample identifier
THREADS="${THREADS:-8}"                    # Number of threads to use

# Tool paths (update these to match your system)
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-/path/to/trimmomatic.jar}"
PICARD_JAR="${PICARD_JAR:-/path/to/picard.jar}"
GATK_JAR="${GATK_JAR:-/path/to/GenomeAnalysisTK.jar}"

# Reference genome and adapter files
REFERENCE_GENOME="${REFERENCE_GENOME:-GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna}" # CpipJ5 assembly (NCBI RefSeq assembly: GCF_015732765.1)
ADAPTER_FILE="${ADAPTER_FILE:-/path/to/adapters/NexteraPE-PE.fa}"

# Input/Output directories
INPUT_DIR="${INPUT_DIR:-./raw_data}"                     # Raw FASTQ files (available from NCBI project PRJNA1209100)
OUTPUT_DIR="${OUTPUT_DIR:-./processed_data}"
QC_DIR="${QC_DIR:-./qc_reports}"
TEMP_DIR="${TEMP_DIR:-./tmp}"

# Create output directories
mkdir -p "${OUTPUT_DIR}" "${QC_DIR}" "${TEMP_DIR}"

#===============================================================================
# STEP 1: INITIAL QUALITY CONTROL
#===============================================================================

# Run FastQC on raw paired-end FASTQ files to assess initial read quality
fastqc "${INPUT_DIR}"/*.fq.gz \
    --threads "${THREADS}" \
    --outdir "${QC_DIR}" \
    --extract

#===============================================================================
# STEP 2: READ TRIMMING AND ADAPTER REMOVAL
#===============================================================================

# Trimmomatic parameters:
# - PE: Paired-end mode
# - ILLUMINACLIP: Remove adapters (seedMismatches:palindromeClipThreshold:simpleClipThreshold)
# - SLIDINGWINDOW: Sliding window trimming (windowSize:requiredQuality)
# - LEADING/TRAILING: Cut bases with quality below threshold from start/end
# - MINLEN: Drop reads shorter than specified length

java -jar "${TRIMMOMATIC_JAR}" PE \
    -threads "${THREADS}" \
    "${INPUT_DIR}/${SAMPLE_NAME}_R1.fastq.gz" \
    "${INPUT_DIR}/${SAMPLE_NAME}_R2.fastq.gz" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}_forward_paired.fq.gz" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}_forward_unpaired.fq.gz" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}_reverse_paired.fq.gz" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}_reverse_unpaired.fq.gz" \
    ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

#===============================================================================
# STEP 3: POST-TRIMMING QUALITY CONTROL
#===============================================================================

# Run FastQC on trimmed files to assess quality improvements
fastqc "${OUTPUT_DIR}"/*_paired.fq.gz \
    --threads "${THREADS}" \
    --outdir "${QC_DIR}" \
    --extract

#===============================================================================
# STEP 4: READ ALIGNMENT
#===============================================================================

# BWA-MEM alignment parameters:
# -t: Number of threads
# -M: Mark shorter split hits as secondary (for Picard compatibility)
# -R: Read group header line

READ_GROUP="@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:${SAMPLE_NAME}_lib"

# Align paired reads to reference genome
bwa mem -t "${THREADS}" -M -R "${READ_GROUP}" \
    "${REFERENCE_GENOME}" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}_forward_paired.fq.gz" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}_reverse_paired.fq.gz" | \
    samtools view -@ "${THREADS}" -bS - | \
    samtools sort -@ "${THREADS}" -o "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam" -

# Index the sorted BAM file
samtools index "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam"

#===============================================================================
# STEP 5: DUPLICATE MARKING
#===============================================================================

# Create directory for duplicate metrics
mkdir -p "${OUTPUT_DIR}/metrics"

# Mark PCR and optical duplicates
java -Xmx8g -jar "${PICARD_JAR}" MarkDuplicates \
    INPUT="${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam" \
    OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.markdup.bam" \
    METRICS_FILE="${OUTPUT_DIR}/metrics/${SAMPLE_NAME}.markdup.metrics.txt" \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true \
    REMOVE_DUPLICATES=false

# Generate mapping statistics
samtools flagstat "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.markdup.bam" \
    > "${OUTPUT_DIR}/metrics/${SAMPLE_NAME}.flagstat.txt"

#===============================================================================
# STEP 6: INDEL REALIGNMENT (GATK3 - Optional for modern pipelines)
#===============================================================================

# Step 6a: Identify target intervals for indel realignment
java -Xmx8g -jar "${GATK_JAR}" \
    -T RealignerTargetCreator \
    -R "${REFERENCE_GENOME}" \
    -I "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.markdup.bam" \
    -o "${TEMP_DIR}/${SAMPLE_NAME}.indel.intervals" \
    -nt 1

# Step 6b: Perform local realignment around indels
java -Xmx8g -jar "${GATK_JAR}" \
    -T IndelRealigner \
    -R "${REFERENCE_GENOME}" \
    -I "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.markdup.bam" \
    -targetIntervals "${TEMP_DIR}/${SAMPLE_NAME}.indel.intervals" \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}.realigned.bam"

# Index the realigned BAM file
samtools index "${OUTPUT_DIR}/${SAMPLE_NAME}.realigned.bam"

#===============================================================================
# STEP 7: COVERAGE CALCULATION
#===============================================================================

# Create output directory for coverage files
mkdir -p "${OUTPUT_DIR}/coverage"

# Calculate coverage statistics from the final processed BAM file
mosdepth \
    --threads "${THREADS}" \
    --by 1000 \
    "${OUTPUT_DIR}/coverage/${SAMPLE_NAME}" \
    "${OUTPUT_DIR}/${SAMPLE_NAME}.realigned.bam"

#===============================================================================
# CLEANUP (OPTIONAL)
#===============================================================================

# Optional: Clean up intermediate files to save disk space
# Uncomment the lines below if you want to remove intermediate files
# rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam"
# rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.markdup.bam"
# rm -rf "${TEMP_DIR}"

#===============================================================================
# PIPELINE COMPLETION SUMMARY
#===============================================================================

# Final processed file: ${OUTPUT_DIR}/${SAMPLE_NAME}.realigned.bam
# Quality reports: ${QC_DIR}/
# Metrics: ${OUTPUT_DIR}/metrics/
# Coverage: ${OUTPUT_DIR}/coverage/