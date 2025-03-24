#!/bin/bash

##############################
### PART 1: Prephase Variants
##############################

### Variables (Update these or pass as arguments)
VARIANTS=$1        # Input filtered biallelic variant file (VCF/BCF)
CHR=$2             # Chromosome or region identifier (e.g., "chr1")
SAMPLE_ID=$3       # Sample identifier
SAMPLE_BAM=$4      # Sample BAM file containing aligned reads

### 1. Extract and Filter Sample-Specific Variants
date
echo "Extracting sample-specific variants for ${SAMPLE_ID} on ${CHR}"

# Extract variants for the sample and filter for genotype quality (GQ) and depth (DP)
bcftools view -Ou -s ${SAMPLE_ID} ${VARIANTS} | \
bcftools filter -S . -Ou -i "FMT/GQ >= 20 & FMT/DP >= 8" \
    -Ou -o /tmp/${SAMPLE_ID}.${CHR}.bcf

# Separate heterozygous calls from non-heterozygous ones
echo "Extracting heterozygous variants..."
bcftools view -Ov -g het /tmp/${SAMPLE_ID}.${CHR}.bcf -o /tmp/${SAMPLE_ID}.${CHR}.hets.vcf
echo "Extracting non-heterozygous variants..."
bcftools view -Ob -g ^het /tmp/${SAMPLE_ID}.${CHR}.bcf -o /tmp/${SAMPLE_ID}.${CHR}.nonhets.bcf

echo "Sample-specific VCF extraction for ${SAMPLE_ID} on ${CHR} completed."
date

### 2. Extract Informative Reads Using extractHAIRS
echo "Extracting informative reads from BAM file..."
extractHAIRS \
    --bam ${SAMPLE_BAM} \
    --VCF /tmp/${SAMPLE_ID}.${CHR}.hets.vcf \
    --out /tmp/${SAMPLE_ID}.${CHR}.hairs
date

### 3. Prephase Heterozygous Variants with HAPCUT2
echo "Running HAPCUT2 for prephasing heterozygous variants..."
HAPCUT2 \
    --fragments /tmp/${SAMPLE_ID}.${CHR}.hairs \
    --vcf /tmp/${SAMPLE_ID}.${CHR}.hets.vcf \
    --out /tmp/${SAMPLE_ID}.${CHR}.hets \
    --outvcf 1
date

### 4. Concatenate Prephased Hets with Non-Hets
echo "Merging prephased heterozygous variants with non-heterozygous variants..."
bcftools concat -Ou /tmp/${SAMPLE_ID}.${CHR}.nonhets.bcf /tmp/${SAMPLE_ID}.${CHR}.hets.phased.VCF | \
bcftools annotate -Ou -x INFO,^FMT/GT,FMT/PS | \
bcftools sort -T /tmp/ -Ob -o /tmp/${SAMPLE_ID}.${CHR}.prephased.bcf

# Index the prephased variant file for downstream use.
bcftools index /tmp/${SAMPLE_ID}.${CHR}.prephased.bcf
date

### 5. Compute Phasing Statistics with Whatshap
# Load Python environment for Whatshap
module load anaconda3/2021.11
source /usr/licensed/anaconda3/2021.11/etc/profile.d/conda.sh
conda activate whatshap

echo "Computing Whatshap statistics..."
mkdir -p whatshap_stats
whatshap stats /tmp/${SAMPLE_ID}.${CHR}.prephased.bcf \
    --tsv /tmp/${SAMPLE_ID}.${CHR}.prephased.whatshap_stats.tsv

# Move stats to a designated folder.
mv /tmp/${SAMPLE_ID}.${CHR}.prephased.whatshap_stats.tsv ./whatshap_stats/

################################################################################
### PART 2: Phasing via SHAPEIT
################################################################################

### Variables for SHAPEIT4 (These can be passed as arguments or defined here)
PREPHASED_VARIANTS_BCF=$1        # Input prephased BCF file (from part 1)
CHR=$2                         # Chromosome/region (same as before)
GENETICMAP=$3                  # Path to the genetic map file
VARIANTS_BCF_HEADER=$(basename ${PREPHASED_VARIANTS_BCF} .prephased.bcf)

# Load Python environment for SHAPEIT4
module load anaconda3/2021.11
source /usr/licensed/anaconda3/2021.11/etc/profile.d/conda.sh
conda activate shapeit4

echo "Running SHAPEIT4 for final phasing..."
shapeit4 \
    --input ${PREPHASED_VARIANTS_BCF} \
    --sequencing \                           # Specify that the data is sequencing-based (not SNP array)
    --use-PS 0.0001 \                        # Expected error rate in phase sets; adjust if needed.
    --effective-size 100000 \                # Effective population size (adjust based on study)
    --region ${CHR} \                        # Target region or chromosome
    --map ${GENETICMAP} \                    # Genetic map file for recombination rates
    --thread 10 \                            # Number of threads for computation
    --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \  # Custom MCMC iterations; modify as needed.
    --pbwt-depth 8 \                         # PBWT depth parameter; adjust if necessary.
    --output /tmp/${VARIANTS_BCF_HEADER}.phased.bcf \  # Output file for phased variants.
    --log ${VARIANTS_BCF_HEADER}.shapeit.log          # Log file for SHAPEIT4 run

# Index the final phased BCF file.
bcftools index /tmp/${VARIANTS_BCF_HEADER}.phased.bcf