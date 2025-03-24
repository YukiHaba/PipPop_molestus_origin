#!/bin/bash

##############################
#### Raw Variant Calling ####
##############################
# Variables used in this section:
#   BAM_LIST: File containing a list of BAM files.
#   REFERENCE_GENOME: Path to the reference genome FASTA file.
#   THREADS: Number of threads for parallel processing.
#   REGION: Genomic region to process (e.g., "chr1:1-20000000").
#   HEADER: Prefix for naming output files.

bcftools mpileup \
    -b ${BAM_LIST} \                                # Use the list of BAM files.
    -f ${REFERENCE_GENOME} \                          # Reference genome file.
    --threads ${THREADS} \                           # Number of threads to use.
    -r ${REGION} \                                  # Target genomic region.
    -Ou -a FORMAT/DP,FORMAT/AD,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR | \
bcftools call \
    --threads ${THREADS} \                           # Use the same thread count.
    -m \                                            # Enable multiallelic calling.
    -a GQ,GP \                                      # Annotate with genotype quality (GQ) and probabilities (GP).
    -Ob -o /tmp/${HEADER}.$(echo ${REGION} | sed 's/:/-/g').allsites.bcf  # Output all-sites BCF file.

# Index the generated BCF file.
bcftools index /tmp/${HEADER}.$(echo ${REGION} | sed 's/:/-/g').allsites.bcf

###################################
#### Variant Quality Filtering ####
###################################
# This section filters raw variant calls to retain high-quality, biallelic sites.
# Criteria for a GOOD site include:
#   1. Being outside repeat regions.
#   2. Not being within 5bp of indels.
#   3. QUAL > 30.
#   4. Mapping Quality (MQ) > 40.
#   5. Less than 25% missing data among samples.
#
# Variables:
#   REGION_FILE: File containing a list of genomic regions.
#   REPEAT_BED: BED file defining repeat regions.
#   HEADER: Prefix for naming files (passed as a command-line argument).
#   SCRIPTS_PATH: Path to the directory containing the filtering worker script.
#   VARIANT_FILE_LIST: A file listing the filtered regional variant files to be concatenated.

# Loop over each region listed in REGION_FILE and process it with the filtering script.
while read -r REGION; do
    bash ${SCRIPTS_PATH}/allsites2good_biallelic.worker.v1.sh \
         ${HEADER}.$(echo ${REGION} | sed 's/:/-/g').allsites.renamed.vcf.gz \  # Input VCF file for this region.
         ${REPEAT_BED} \                   # Provide the repeat regions file.
         ${HEADER}.$(echo ${REGION} | sed 's/:/-/g') &   # Output file prefix.
done < ${REGION_FILE} # 20Mb chunks

# Concatenate the regional filtered VCF files into a single file.
bcftools concat \
    --allow-overlaps \
    --file-list ${VARIANT_FILE_LIST} \   # Text file listing all regional VCF files.
    --threads 20 \                       # Use 20 threads for the operation.
    -Oz -o /tmp/${HEADER}_combined.vcf.gz  # Output concatenated VCF file.
    
# Index the concatenated VCF.
bcftools index /tmp/${HEADER}_combined.vcf.gz

# Extract variants that fall in accessible genomic regions.
# Variables:
#   CHR: Chromosome identifier.
#   ACCESSIBLE_BED: BED file with accessible regions for the chromosome.
#   INPUT_BCF: Name of the input BCF file containing filtered variants.
CHR=<chromosome_ID>                              # Replace with your chromosome identifier.
ACCESSIBLE_BED="/path/to/accessible_regions/chr${CHR}.bed"  # Update path to accessible regions BED file.
INPUT_BCF="${HEADER}.chr${CHR}.filtered.bcf"     # Adjust input file naming as needed.

bcftools view \
    ${INPUT_BCF} \                       # Input BCF file with filtered variants.
    -R ${ACCESSIBLE_BED} \               # Restrict to accessible regions.
    -Oz -o ${HEADER}.chr${CHR}.accessible.vcf.gz &  # Output VCF with accessible variants.