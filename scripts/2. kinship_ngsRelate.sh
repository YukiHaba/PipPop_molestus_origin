#!/bin/bash

################################################################################
# SCRIPT: 2. Kinship Coefficient Calculation
# AUTHOR: Yuki Haba
# REFERENCE: Yuki Haba et al. 2025 Ancient origin of an urban underground 
#            mosquito (PMID: 39975080)
# DESCRIPTION: Calculate kinship coefficients and relatedness estimates between
#              individuals using ngsRelate from high-quality variant data.
################################################################################

# Exit on any error
set -e

#===============================================================================
# CONFIGURATION VARIABLES - MODIFY THESE PATHS FOR YOUR SYSTEM
#===============================================================================

# Input files
INPUT_VCF="${INPUT_VCF:-good_biallelic_snps.rm.combined.accessible.vcf.gz}" # High-quality variant file from script 1
INDIVIDUAL_LIST="${INDIVIDUAL_LIST:-/path/to/individuals.txt}"               # File with individual names

# Analysis parameters
THREADS="${THREADS:-20}"                                     # Number of processing threads
N_INDIVIDUALS="${N_INDIVIDUALS:-840}"                        # Total number of individuals

# Output settings
OUTPUT_PREFIX="${OUTPUT_PREFIX:-ngsRelate_results}"          # Output file prefix
OUTPUT_DIR="${OUTPUT_DIR:-./kinship_analysis}"              

# Create output directory
mkdir -p "${OUTPUT_DIR}"

#===============================================================================
# STEP 1: KINSHIP COEFFICIENT CALCULATION
#===============================================================================

# Calculate pairwise kinship coefficients using ngsRelate
# This analysis estimates relatedness between all pairs of individuals
# using genome-wide SNP data

ngsRelate \
    --haps "${INPUT_VCF}" \
    --names "${INDIVIDUAL_LIST}" \
    --nthreads "${THREADS}" \
    --nind "${N_INDIVIDUALS}" \
    --outprefix "${OUTPUT_DIR}/${OUTPUT_PREFIX}"

#===============================================================================
# PIPELINE COMPLETION SUMMARY
#===============================================================================

# Output files:
# - Kinship results: ${OUTPUT_DIR}/${OUTPUT_PREFIX}