#!/bin/bash

################################################################################
# SCRIPT: 1. Variant Calling and Filtering Pipeline
# AUTHOR: Yuki Haba
# REFERENCE: Yuki Haba et al. 2025 Ancient origin of an urban underground 
#            mosquito (PMID: 39975080)
# DESCRIPTION: Pipeline for variant calling from aligned BAM files,
#              quality filtering, repeat masking, accessible site extraction,
#              and linkage disequilibrium pruning for population genomics.
################################################################################

# Exit on any error
set -e

#===============================================================================
# CONFIGURATION VARIABLES - MODIFY THESE PATHS FOR YOUR SYSTEM
#===============================================================================

# Input files and directories
BAM_LIST="${BAM_LIST:-/path/to/bam_list.txt}"                    # File containing list of realigned BAM files from script 0
REFERENCE_GENOME="${REFERENCE_GENOME:-GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna}" # CpipJ5 assembly (NCBI RefSeq assembly: GCF_015732765.1)
REPEAT_BED="${REPEAT_BED:-GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_rm.out.gz}" # Repeat regions file (provided in repository)
ACCESSIBLE_BED="${ACCESSIBLE_BED:-/path/to/accessible_sites.bed}" # Accessible genomic regions BED file

# Processing parameters
THREADS="${THREADS:-20}"                                         # Number of threads for parallel processing
REGION="${REGION:-chr1:1-20000000}"                             # Genomic region to process (20Mb chunks recommended)
HEADER="${HEADER:-variant_calling_run}"                         # Prefix for naming output files

# Output directories
OUTPUT_DIR="${OUTPUT_DIR:-./variant_calling_output}"
TEMP_DIR="${TEMP_DIR:-./tmp}"

# Create output directories
mkdir -p "${OUTPUT_DIR}" "${TEMP_DIR}"

#===============================================================================
# STEP 1: RAW VARIANT CALLING
#===============================================================================

# Call variants using bcftools mpileup and call
# Process in chunks for memory efficiency and parallelization
bcftools mpileup \
    --bam-list "${BAM_LIST}" \
    --fasta-ref "${REFERENCE_GENOME}" \
    --threads "${THREADS}" \
    --regions "${REGION}" \
    --output-type u \
    --annotate FORMAT/DP,FORMAT/AD,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR | \
bcftools call \
    --threads "${THREADS}" \
    --multiallelic-caller \
    --annotate GQ,GP \
    --output-type b \
    --output "${TEMP_DIR}/${HEADER}.$(echo ${REGION} | sed 's/:/-/g').allsites.bcf"

# Index the generated BCF file
bcftools index "${TEMP_DIR}/${HEADER}.$(echo ${REGION} | sed 's/:/-/g').allsites.bcf"

#===============================================================================
# STEP 2: VARIANT QUALITY FILTERING
#===============================================================================

# Filter variants to retain high-quality, biallelic sites
# Quality criteria for a GOOD site:
#   1. Outside repeat regions (using repeat mask BED file)
#   2. Not within 5bp of indels
#   3. QUAL > 30 (variant quality score)
#   4. Mapping Quality (MQ) > 40
#   5. Less than 25% missing data among samples
#   6. Biallelic SNPs only

# Set input file variable (this should match output from Step 1)
ALLSITES_VARIANTS="${TEMP_DIR}/${HEADER}.$(echo ${REGION} | sed 's/:/-/g').allsites.bcf"

# Apply filtering pipeline
bcftools filter \
    --SnpGap 5 \
    "${ALLSITES_VARIANTS}" \
    --output-type u | \
bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    --types snps \
    --output-type u | \
bcftools view \
    --output-type b \
    --targets-file ^"${REPEAT_BED}" \
    --include 'F_MISSING < 0.25 & MQ > 40 & QUAL > 30' \
    --output "${TEMP_DIR}/${REGION}.good_biallelic_snps.rm.bcf"

# Index the filtered BCF file
bcftools index "${TEMP_DIR}/${REGION}.good_biallelic_snps.rm.bcf"

#===============================================================================
# STEP 3: CONCATENATE REGIONAL VCF FILES
#===============================================================================

# Concatenate filtered VCF files from different regions into a single file
# VARIANT_FILE_LIST should contain paths to all regional filtered BCF files
VARIANT_FILE_LIST="${VARIANT_FILE_LIST:-/path/to/variant_file_list.txt}"

bcftools concat \
    --allow-overlaps \
    --file-list "${VARIANT_FILE_LIST}" \
    --threads "${THREADS}" \
    --output-type z \
    --output "${TEMP_DIR}/good_biallelic_snps.rm.combined.vcf.gz"

# Index the concatenated VCF
bcftools index "${TEMP_DIR}/good_biallelic_snps.rm.combined.vcf.gz"

#===============================================================================
# STEP 4: EXTRACT ACCESSIBLE SITES
#===============================================================================

# Extract variants that fall within accessible genomic regions
INPUT_VCF="${TEMP_DIR}/good_biallelic_snps.rm.combined.vcf.gz"

bcftools view \
    "${INPUT_VCF}" \
    --regions-file "${ACCESSIBLE_BED}" \
    --output-type z \
    --output "${OUTPUT_DIR}/good_biallelic_snps.rm.combined.accessible.vcf.gz"

# Index the accessible sites VCF
bcftools index "${OUTPUT_DIR}/good_biallelic_snps.rm.combined.accessible.vcf.gz"


#===============================================================================
# STEP 5: LINKAGE DISEQUILIBRIUM PRUNING (FOR DOWNSTREAM ANALYSES)
#===============================================================================

# Generate unlinked SNPs for analyses that require independent markers
# (e.g., PCA, ADMIXTURE, some phylogenetic analyses)

VCFHEADER="good_biallelic_snps.rm.combined.accessible"
INPUT_VCF="${OUTPUT_DIR}/${VCFHEADER}.vcf.gz"

# Step 5a: Convert VCF to PLINK map/ped format
plink --vcf "${INPUT_VCF}" \
    --allow-extra-chr \
    --recode \
    --out "${TEMP_DIR}/${VCFHEADER}" \
    --set-missing-var-ids @:#[Cpip] \
    --threads "${THREADS}"

# Step 5b: Identify SNPs in linkage disequilibrium for removal
# Parameters: window size (200 SNPs), step size (20 SNPs), r² threshold (0.2)
# Remove SNPs with r² > 0.2 within sliding windows
plink --file "${TEMP_DIR}/${VCFHEADER}" \
    --allow-extra-chr \
    --indep-pairwise 200 20 0.2 \
    --out "${TEMP_DIR}/${VCFHEADER}" \
    --threads "${THREADS}"

# Step 5c: Extract pruned SNPs and convert back to VCF format
plink --file "${TEMP_DIR}/${VCFHEADER}" \
    --allow-extra-chr \
    --extract "${TEMP_DIR}/${VCFHEADER}.prune.in" \
    --out "${TEMP_DIR}/${VCFHEADER}.LDpruned" \
    --recode vcf \
    --threads "${THREADS}"

# Step 5d: Compress, rename samples (if needed), and index final VCF
bgzip "${TEMP_DIR}/${VCFHEADER}.LDpruned.vcf"

# Optional: Rename samples in the VCF file if sample name file exists
SAMPLE_NAMES_FILE="${SAMPLE_NAMES_FILE:-/path/to/sample_names.txt}"
if [[ -f "${SAMPLE_NAMES_FILE}" ]]; then
    bcftools reheader \
        --samples "${SAMPLE_NAMES_FILE}" \
        "${TEMP_DIR}/${VCFHEADER}.LDpruned.vcf.gz" \
        --output "${OUTPUT_DIR}/${VCFHEADER}.LDpruned.vcf.gz"
else
    mv "${TEMP_DIR}/${VCFHEADER}.LDpruned.vcf.gz" "${OUTPUT_DIR}/${VCFHEADER}.LDpruned.vcf.gz"
fi

# Index the final pruned VCF
bcftools index "${OUTPUT_DIR}/${VCFHEADER}.LDpruned.vcf.gz"

#===============================================================================
# CLEANUP (OPTIONAL)
#===============================================================================

# Optional: Clean up intermediate files to save disk space
# Uncomment the lines below if you want to remove intermediate files
# rm -f "${TEMP_DIR}"/*.bcf
# rm -f "${TEMP_DIR}"/*.vcf.gz
# rm -f "${TEMP_DIR}"/*.ped
# rm -f "${TEMP_DIR}"/*.map

#===============================================================================
# PIPELINE COMPLETION SUMMARY
#===============================================================================

# Final output files:
# - All accessible variants: ${OUTPUT_DIR}/good_biallelic_snps.rm.combined.accessible.vcf.gz
# - LD-pruned variants: ${OUTPUT_DIR}/${VCFHEADER}.LDpruned.vcf.gz
# - Intermediate files: ${TEMP_DIR}/
