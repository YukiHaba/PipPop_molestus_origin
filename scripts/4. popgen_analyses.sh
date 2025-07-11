#!/bin/bash

################################################################################
# SCRIPT: 4. Population Genomics Analysis Pipeline
# AUTHOR: Yuki Haba
# REFERENCE: Yuki Haba et al. 2025 Ancient origin of an urban underground 
#            mosquito (PMID: 39975080)
# DESCRIPTION: Population genetic analyses including PCA, f-statistics,
#              nucleotide diversity, and demographic history inference 
################################################################################

# Exit on any error
set -e

#===============================================================================
# CONFIGURATION VARIABLES - MODIFY THESE PATHS FOR YOUR SYSTEM
#===============================================================================

# Input files
BEAGLE_FILE="${BEAGLE_FILE:-good_biallelic_snps.rm.combined.accessible.LDpruned.beagle.gz}" # Beagle format from script 1
VCF_FILE="${VCF_FILE:-good_biallelic_snps.rm.combined.accessible.vcf.gz}"                   # VCF from script 1
LDPRUNED_VCF="${LDPRUNED_VCF:-good_biallelic_snps.rm.combined.accessible.LDpruned.vcf.gz}" # LD-pruned from script 1
ALLSITES_BCF="${ALLSITES_BCF:-/path/to/allsites.accessible.bcf}"                            # All sites BCF with invariants
POPULATION_LIST="${POPULATION_LIST:-/path/to/populations.list}"                             # Population assignment file
POPULATION_TREE="${POPULATION_TREE:-/path/to/tree.nwk}"                                    # Population tree in Newick format

# Analysis parameters
THREADS="${THREADS:-20}"                                      # Number of threads
HEADER="${HEADER:-popgen_analysis}"                          # Output file prefix
MIN_MAF="${MIN_MAF:-0.05}"                                   # Minimum minor allele frequency

# Tool environments
PCANGSD_ENV="${PCANGSD_ENV:-pcangsd}"                        # Conda environment for PCAngsd
TREEMIX_ENV="${TREEMIX_ENV:-treemix}"                        # Conda environment for Treemix
PIXY_ENV="${PIXY_ENV:-pixy}"                                 # Conda environment for pixy

# Output directories
OUTPUT_DIR="${OUTPUT_DIR:-./popgen_results}"
TEMP_DIR="${TEMP_DIR:-./tmp}"

# Create output directories
mkdir -p "${OUTPUT_DIR}" "${TEMP_DIR}"

#===============================================================================
# STEP 1: PRINCIPAL COMPONENT ANALYSIS (PCA)
#===============================================================================

# Activate PCAngsd environment
conda activate "${PCANGSD_ENV}"

# Run PCA with SNP loadings and site information
pcangsd \
    --beagle "${BEAGLE_FILE}" \
    --threads "${THREADS}" \
    --minMaf "${MIN_MAF}" \
    --snp_weights \
    --sites_save \
    --out "${OUTPUT_DIR}/${HEADER}_pca"

conda deactivate

#===============================================================================
# STEP 2: F3 STATISTICS (THREE-POPULATION TEST)
#===============================================================================

# Activate Treemix environment
conda activate "${TREEMIX_ENV}"

# Convert VCF to Treemix format
populations \
    --in-vcf "${LDPRUNED_VCF}" \
    --popmap "${POPULATION_LIST}" \
    --threads "${THREADS}" \
    --treemix \
    --out-path "${OUTPUT_DIR}"

# Run three-population test with block size of 500 SNPs
threepop \
    --infile "${OUTPUT_DIR}/good_biallelic_snps.rm.combined.accessible.LDpruned.treemix.gz" \
    --block-size 500 \
    > "${OUTPUT_DIR}/${HEADER}.threepop.out"

conda deactivate

#===============================================================================
# STEP 3: F4 STATISTICS AND D-STATISTICS
#===============================================================================

# Run D-statistics analysis using DtriosParallel
DtriosParallel \
    --keep-intermediate \
    --tree "${POPULATION_TREE}" \
    --cores "${THREADS}" \
    "${POPULATION_LIST}" \
    "${LDPRUNED_VCF}"

# Clean up temporary file naming
rename '_tmp.popset_file' '' "${TEMP_DIR}"/DTparallel_*

# Run Fbranch analysis for tree-based statistics
Dsuite Fbranch \
    "${POPULATION_TREE}" \
    "${TEMP_DIR}/DTparallel_combined_tree.txt" \
    > "${OUTPUT_DIR}/Fbranch.DTparallel_combined_tree.txt"

# Process Fbranch results for visualization
DSUITE_UTILS="${DSUITE_UTILS:-/path/to/Dsuite/utils}"
python3 "${DSUITE_UTILS}/dtools.py" \
    "${OUTPUT_DIR}/Fbranch.DTparallel_combined_tree.txt" \
    "${POPULATION_TREE}" \
    --ladderize \
    --color-cutoff 0.5 \
    --tree-label-size 3 \
    --run-name "${OUTPUT_DIR}/Fbranch.DTparallel"

#===============================================================================
# STEP 4: NUCLEOTIDE DIVERSITY ANALYSIS
#===============================================================================

# Activate pixy environment
conda activate "${PIXY_ENV}"

# Calculate nucleotide diversity (pi) and divergence (dxy) statistics
# Note: Requires all-sites BCF file containing invariant sites
pixy --stats dxy \
    --vcf "${ALLSITES_BCF}" \
    --populations "${POPULATION_LIST}" \
    --n_cores "${THREADS}" \
    --output_folder "${OUTPUT_DIR}"

conda deactivate

#===============================================================================
# STEP 5: DEMOGRAPHIC HISTORY INFERENCE (MSMC2)
#===============================================================================

#===============================================================================
# STEP 5A: MSMC2 PREPARATION
#===============================================================================

# MSMC2 analysis variables (modify as needed)
SAMPLE_ID="${SAMPLE_ID:-sample_name}"                                       # Sample identifier
CHR="${CHR:-chr1}"                                                          # Chromosome identifier
SAMPLE_BAM="${SAMPLE_BAM:-/path/to/sample.realigned.bam}"                  # Sample BAM file from script 0
PHASED_VCF="${PHASED_VCF:-${SAMPLE_ID}.${CHR}.phased.bcf}"                 # Phased variant file from script 3
ACCESSIBLE_BED="${ACCESSIBLE_BED:-/path/to/accessible.bed}"                 # Accessible regions BED
REFERENCE="${REFERENCE:-GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna}" # CpipJ5 assembly (NCBI RefSeq assembly: GCF_015732765.1)

# Calculate mean sequencing depth for the sample
DEPTH=$(samtools depth -b "${ACCESSIBLE_BED}" "${SAMPLE_BAM}" | awk '{sum += $3} END {print sum / NR}')

# Create sample-specific phased VCF
bcftools view \
    --output-type z \
    --samples "${SAMPLE_ID}" \
    "${PHASED_VCF}" \
    > "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.phased.vcf.gz"

# Generate mask file using mpileup and bamCaller.py
BAMCALLER="${BAMCALLER:-/path/to/msmc-tools/bamCaller.py}"
samtools mpileup \
    -B -q20 -Q20 -C50 -g \
    -l "${ACCESSIBLE_BED}" \
    -f "${REFERENCE}" \
    "${SAMPLE_BAM}" | \
bcftools call -c -V indels | \
"${BAMCALLER}" "${DEPTH}" "${TEMP_DIR}/${SAMPLE_ID}.mask.chr${CHR}.bed.gz" | \
bgzip -c > "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.vcf.gz"

# Index VCF files
bcftools index "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.vcf.gz"
bcftools index "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.phased.vcf.gz"

# Merge phased variants with mask variants
bcftools merge \
    --force-samples \
    "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.vcf.gz" \
    "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.phased.vcf.gz" | \
awk 'BEGIN {OFS="\t"}
     $0 ~ /^##/ {print $0}
     $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
     $0 !~ /^#/ {if(substr($11, 1, 3) != "./.") $10 = $11; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' | \
bcftools view -Oz > "${TEMP_DIR}/${SAMPLE_ID}.chr${CHR}.final.vcf.gz"

#===============================================================================
# STEP 5B: MSMC2 MAIN ANALYSIS
#===============================================================================

# MSMC2 analysis variables
RUN_NAME="${RUN_NAME:-msmc_analysis}"                        # Analysis run identifier
SAMPLE1="${SAMPLE1:-sample1}"                                # Individual sample IDs
SAMPLE2="${SAMPLE2:-sample2}"
SAMPLE3="${SAMPLE3:-sample3}"
SAMPLE4="${SAMPLE4:-sample4}"
POP1="${POP1:-population1}"                                  # Population identifiers
POP2="${POP2:-population2}"
MSMC_TOOLS="${MSMC_TOOLS:-/path/to/msmc-tools}"             # Path to MSMC tools

# Create MSMC output directory
mkdir -p "${OUTPUT_DIR}/${RUN_NAME}"

# Generate multihetsep files for each chromosome
for CHR in 1 2 3; do
    # Mask files for each sample
    MASKFILE_SAMPLE1="${SAMPLE1}.mask.chr${CHR}.bed.gz"
    MASKFILE_SAMPLE2="${SAMPLE2}.mask.chr${CHR}.bed.gz"
    MASKFILE_SAMPLE3="${SAMPLE3}.mask.chr${CHR}.bed.gz"
    MASKFILE_SAMPLE4="${SAMPLE4}.mask.chr${CHR}.bed.gz"
    
    # Final VCF files from MSMC prep
    FINALVCF1="${SAMPLE1}.chr${CHR}.final.vcf.gz"
    FINALVCF2="${SAMPLE2}.chr${CHR}.final.vcf.gz"
    FINALVCF3="${SAMPLE3}.chr${CHR}.final.vcf.gz"
    FINALVCF4="${SAMPLE4}.chr${CHR}.final.vcf.gz"
    
    # Generate multihetsep file
    python3 "${MSMC_TOOLS}/generate_multihetsep.py" \
        --chr "chr${CHR}" \
        --mask "${ACCESSIBLE_BED}" \
        --mask "${TEMP_DIR}/${MASKFILE_SAMPLE1}" \
        --mask "${TEMP_DIR}/${MASKFILE_SAMPLE2}" \
        --mask "${TEMP_DIR}/${MASKFILE_SAMPLE3}" \
        --mask "${TEMP_DIR}/${MASKFILE_SAMPLE4}" \
        "${TEMP_DIR}/${FINALVCF1}" \
        "${TEMP_DIR}/${FINALVCF2}" \
        "${TEMP_DIR}/${FINALVCF3}" \
        "${TEMP_DIR}/${FINALVCF4}" \
        > "${OUTPUT_DIR}/${RUN_NAME}/${RUN_NAME}.chr${CHR}.multihetsep.txt" &
done

wait

# Run MSMC2 analyses
# Within-population coalescent analyses
msmc2 -t 4 -s -I 0-1,0-2,0-3,1-2,1-3,2-3 \
    -o "${OUTPUT_DIR}/${RUN_NAME}/${POP1}" \
    "${OUTPUT_DIR}/${RUN_NAME}"/*multihetsep*.txt &

msmc2 -t 4 -s -I 4-5,4-6,4-7,5-6,5-7,6-7 \
    -o "${OUTPUT_DIR}/${RUN_NAME}/${POP2}" \
    "${OUTPUT_DIR}/${RUN_NAME}"/*multihetsep*.txt &

# Cross-population coalescent analysis
msmc2 -t 8 -s -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 \
    -o "${OUTPUT_DIR}/${RUN_NAME}/cross" \
    "${OUTPUT_DIR}/${RUN_NAME}"/*multihetsep*.txt &

wait

# Combine cross-coalescence results
python3 "${MSMC_TOOLS}/combineCrossCoal.py" \
    "${OUTPUT_DIR}/${RUN_NAME}/cross.final.txt" \
    "${OUTPUT_DIR}/${RUN_NAME}/${POP1}.final.txt" \
    "${OUTPUT_DIR}/${RUN_NAME}/${POP2}.final.txt" \
    > "${OUTPUT_DIR}/${RUN_NAME}/${RUN_NAME}.combined.final.txt"

#===============================================================================
# STEP 5C: MSMC2 BOOTSTRAPPING (OPTIONAL)
#===============================================================================

# Uncomment the following section to perform bootstrapping for confidence intervals
# N_BOOT="${N_BOOT:-100}"
# python3 "${MSMC_TOOLS}/multihetsep_bootstrap.py" \
#     --n_bootstrap "${N_BOOT}" \
#     --seed 10000000 \
#     --chunks_per_chromosome 20 \
#     --nr_chromosomes 3 \
#     "${OUTPUT_DIR}/${RUN_NAME}/boot" \
#     "${OUTPUT_DIR}/${RUN_NAME}/${RUN_NAME}.chr1.multihetsep.txt" \
#     "${OUTPUT_DIR}/${RUN_NAME}/${RUN_NAME}.chr2.multihetsep.txt" \
#     "${OUTPUT_DIR}/${RUN_NAME}/${RUN_NAME}.chr3.multihetsep.txt"

#===============================================================================
# CLEANUP (OPTIONAL)
#===============================================================================

# Optional: Clean up intermediate files to save disk space
# Uncomment the lines below if you want to remove intermediate files
# rm -f "${TEMP_DIR}"/*.vcf.gz
# rm -f "${TEMP_DIR}"/*.bed.gz

#===============================================================================
# PIPELINE COMPLETION SUMMARY
#===============================================================================

# Final output files:
# - PCA results: ${OUTPUT_DIR}/${HEADER}_pca.*
# - f3 statistics: ${OUTPUT_DIR}/${HEADER}.threepop.out
# - D-statistics: ${OUTPUT_DIR}/Fbranch.DTparallel_combined_tree.txt
# - Nucleotide diversity: ${OUTPUT_DIR}/pixy_*
# - MSMC2 results: ${OUTPUT_DIR}/${RUN_NAME}/