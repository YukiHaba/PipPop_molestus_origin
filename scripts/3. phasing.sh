#!/bin/bash

################################################################################
# SCRIPT: 3. Haplotype Phasing Pipeline
# AUTHOR: Yuki Haba
# REFERENCE: Yuki Haba et al. 2025 Ancient origin of an urban underground 
#            mosquito (PMID: 39975080)
# DESCRIPTION: Two-stage haplotype phasing pipeline using HAPCUT2 for initial
#              prephasing followed by statistical phasing. 
#              Includes quality filtering and phasing statistics.
################################################################################

# Exit on any error
set -e

#===============================================================================
# CONFIGURATION VARIABLES - MODIFY THESE PATHS FOR YOUR SYSTEM
#===============================================================================

# Input files
VARIANTS="${VARIANTS:-good_biallelic_snps.rm.combined.accessible.vcf.gz}" # Filtered biallelic variant file from script 1
SAMPLE_ID="${SAMPLE_ID:-sample_name}"                                     # Sample identifier
SAMPLE_BAM="${SAMPLE_BAM:-/path/to/sample.bam}"                          # Sample BAM file from script 0
CHR="${CHR:-chr1}"                                                        # Chromosome/region identifier

# Phasing parameters
GENETIC_MAP="${GENETIC_MAP:-/path/to/genetic_map.txt}"       # Genetic map file for recombination rates

# Tool paths and environments
WHATSHAP_ENV="${WHATSHAP_ENV:-whatshap}"                     # Conda environment for whatshap
SHAPEIT4_ENV="${SHAPEIT4_ENV:-shapeit4}"                     # Conda environment for SHAPEIT4

# Output directories
OUTPUT_DIR="${OUTPUT_DIR:-./phasing_output}"
TEMP_DIR="${TEMP_DIR:-./tmp}"
STATS_DIR="${STATS_DIR:-./phasing_stats}"

# Create output directories
mkdir -p "${OUTPUT_DIR}" "${TEMP_DIR}" "${STATS_DIR}"

#===============================================================================
# PART 1: VARIANT PREPHASING WITH HAPCUT2
#===============================================================================

#===============================================================================
# STEP 1: EXTRACT AND FILTER SAMPLE-SPECIFIC VARIANTS
#===============================================================================

# Extract variants for the specific sample and apply quality filters
# Filters: Genotype Quality (GQ) >= 20, Depth (DP) >= 8
bcftools view \
    --output-type u \
    --samples "${SAMPLE_ID}" \
    "${VARIANTS}" | \
bcftools filter \
    --set-GTs . \
    --output-type u \
    --include "FMT/GQ >= 20 & FMT/DP >= 8" \
    --output "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.bcf"

# Separate heterozygous from non-heterozygous variants
bcftools view \
    --output-type v \
    --genotype het \
    "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.bcf" \
    --output "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hets.vcf"

bcftools view \
    --output-type b \
    --genotype ^het \
    "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.bcf" \
    --output "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.nonhets.bcf"

#===============================================================================
# STEP 2: EXTRACT INFORMATIVE READS
#===============================================================================

# Extract reads that span multiple heterozygous sites for phasing
extractHAIRS \
    --bam "${SAMPLE_BAM}" \
    --VCF "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hets.vcf" \
    --out "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hairs"

#===============================================================================
# STEP 3: PREPHASE HETEROZYGOUS VARIANTS WITH HAPCUT2
#===============================================================================

# Phase heterozygous variants using read-based information
HAPCUT2 \
    --fragments "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hairs" \
    --vcf "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hets.vcf" \
    --out "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hets" \
    --outvcf 1

#===============================================================================
# STEP 4: MERGE PREPHASED AND NON-HETEROZYGOUS VARIANTS
#===============================================================================

# Combine prephased heterozygous variants with non-heterozygous variants
bcftools concat \
    --output-type u \
    "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.nonhets.bcf" \
    "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hets.phased.VCF" | \
bcftools annotate \
    --output-type u \
    --remove INFO,^FMT/GT,FMT/PS | \
bcftools sort \
    --temp-dir "${TEMP_DIR}" \
    --output-type b \
    --output "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.prephased.bcf"

# Index the prephased variant file
bcftools index "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.prephased.bcf"

#===============================================================================
# STEP 5: CALCULATE PHASING STATISTICS
#===============================================================================

# Activate whatshap environment and compute phasing statistics
conda activate "${WHATSHAP_ENV}"

whatshap stats \
    "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.prephased.bcf" \
    --tsv "${STATS_DIR}/${SAMPLE_ID}.${CHR}.prephased.whatshap_stats.tsv"

conda deactivate

#===============================================================================
# PART 2: STATISTICAL PHASING WITH SHAPEIT4
#===============================================================================

#===============================================================================
# STEP 6: FINAL PHASING WITH SHAPEIT4
#===============================================================================

# Activate SHAPEIT4 environment
conda activate "${SHAPEIT4_ENV}"

PREPHASED_BCF="${TEMP_DIR}/${SAMPLE_ID}.${CHR}.prephased.bcf"
VARIANTS_BCF_HEADER=$(basename "${PREPHASED_BCF}" .prephased.bcf)

# Run SHAPEIT4 for final statistical phasing using genetic maps
shapeit4 \
    --input "${PREPHASED_BCF}" \
    --sequencing \
    --use-PS 0.0001 \
    --effective-size 100000 \
    --region "${CHR}" \
    --map "${GENETIC_MAP}" \
    --thread 10 \
    --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
    --pbwt-depth 8 \
    --output "${OUTPUT_DIR}/${VARIANTS_BCF_HEADER}.phased.bcf" \
    --log "${OUTPUT_DIR}/${VARIANTS_BCF_HEADER}.shapeit.log"

# Index the final phased BCF file
bcftools index "${OUTPUT_DIR}/${VARIANTS_BCF_HEADER}.phased.bcf"

conda deactivate

#===============================================================================
# CLEANUP (OPTIONAL)
#===============================================================================

# Optional: Clean up intermediate files to save disk space
# Uncomment the lines below if you want to remove intermediate files
# rm -f "${TEMP_DIR}/${SAMPLE_ID}.${CHR}".{bcf,hets.vcf,nonhets.bcf,hairs}
# rm -f "${TEMP_DIR}/${SAMPLE_ID}.${CHR}.hets.phased.VCF"

#===============================================================================
# PIPELINE COMPLETION SUMMARY
#===============================================================================

# Final output files:
# - Prephased variants: ${TEMP_DIR}/${SAMPLE_ID}.${CHR}.prephased.bcf
# - Final phased variants: ${OUTPUT_DIR}/${SAMPLE_ID}.${CHR}.phased.bcf
# - Phasing statistics: ${STATS_DIR}/${SAMPLE_ID}.${CHR}.prephased.whatshap_stats.tsv
# - SHAPEIT4 log: ${OUTPUT_DIR}/${SAMPLE_ID}.${CHR}.shapeit.log