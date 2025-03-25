#!/bin/bash

##############################
# PCA
##############################
# Load Python environment for PCANGSD
conda activate pcangsd

# Run PCAngsd with SNP loadings and site information saved
pcangsd \
    --beagle good_biallelic_snps.rm.bcf_combined.accessible.LDpruned.beagle.gz \
    --threads 5 \
    --minMaf 0.05 \
    --snp_weights \
    --sites_save \
    -o ${HEADER}

##############################
# f3
##############################
# Load Python environment for Treemix
conda activate treemix

HEADER=good_biallelic_snps.rm.bcf_combined.accessible.LDpruned

# get treemix file
populations -V $HEADER.vcf.gz \
-M populations.list \
-t 20 --treemix --no-hap-exports \
-O .

# Run three-population test with a block size of 500 SNPs
threepop -i $HEADER.treemix.gz \
    -k 500 \
    > ${HEADER}.threepop.out

##############################
# f4/D and Fbranch
##############################
# Run DtriosParallel with the intermediate files preserved
DtriosParallel \
    --keep-intermediate \
    -t population_tree.nwk \
    --cores 20 \
    population.list \
    good_biallelic_snps.rm.bcf_combined.accessible.LDpruned.vcf.gz

# Rename temporary output files (remove unwanted string from filenames)
rename _tmp.popset_file "" /tmp/DTparallel_*

# Run Fbranch analysis using Dsuite with the combined tree file
Dsuite Fbranch \
    population_tree.nwk \
    /tmp/DTparallel_combined_tree.txt \
    > /tmp/Fbranch.DTparallel_combined_tree.txt

# Process Fbranch results for visualization
python3 /path/to/Dsuite/utils/dtools.py \
    Fbranch.DTparallel_combined_tree.txt \
    population_tree.nwk \
    --ladderize \
    --color-cutoff 0.5 \
    --tree-label-size 3 \
    --run-name Fbranch.DTparallel

##############################
# Nucleotide Diversity
##############################
# Load Python environment for pixy
conda activate pixy

# Run pixy to calculate dxy statistics. NOTE: use allsite bcf that contains invariant sites
pixy --stats dxy \
    --vcf PipPop.allsites.accessible.bcf \
    --populations population.list \
    --n_cores 20 \
    --output_folder . 

################################################################################
# MSMC 
################################################################################
### MSMC Prep
# Variables for MSMC prep
$SAMPLE                          # Sample identifier
$CHR                             # Chromosome identifier (e.g., "NC_051861.1")
BAM=/path/to/$SAMPLE.bam               # BAM file for the sample
PHASED_GENOME=/path/to/good_biallelic_snps.rm.bcf_combined.accessible.phased.bcf  # Input phased VCF file
REGION_TO_INCLUDE_perCHR_BED=/path/to/PipPop.accessible_sites.bed  # Accessible regions BED file
REF=/path/to/reference.fasta           # Reference genome file

# Calculate mean depth over the included regions
DEPTH=$(samtools depth -b ${REGION_TO_INCLUDE_perCHR_BED} ${BAM} | awk '{sum += $3} END {print sum / NR}')

# Create a phased VCF file for a given sample and chromosome
bcftools view -Oz -s ${SAMPLE} ${PHASED_GENOME} \
    > /tmp/${SAMPLE}.chr${CHR}.phased.vcf.gz

# Generate a sample-specific mask file using mpileup, variant calling, and bamCaller.py
samtools mpileup -B -q20 -Q20 -C50 -g \
    -l ${REGION_TO_INCLUDE_perCHR_BED} -f ${REF} ${BAM} | \
bcftools call -c -V indels | \
/path/to/msmc-tools/bamCaller.py ${DEPTH} /tmp/${SAMPLE}.mask.chr${CHR}.bed.gz | \
bgzip -c > /tmp/${SAMPLE}.chr${CHR}.vcf.gz

# Index the VCF files
bcftools index /tmp/${SAMPLE}.chr${CHR}.vcf.gz
bcftools index /tmp/${SAMPLE}.chr${CHR}.phased.vcf.gz

# Merge phased variants with mask variants into a final VCF for MSMC input
bcftools merge --force-samples /tmp/${SAMPLE}.chr${CHR}.vcf.gz /tmp/${SAMPLE}.chr${CHR}.phased.vcf.gz | \
awk 'BEGIN {OFS="\t"}
     $0 ~ /^##/ {print $0}
     $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
     $0 !~ /^#/ {if(substr($11, 1, 3) != "./.") $10 = $11; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' | \
bcftools view -Oz > /tmp/${SAMPLE}.chr${CHR}.final.vcf.gz

### Main MSMC2 Run
# Variables for MSMC2 main run (and bootstrapping)
RUN_NAME=msmc_run_example         # Unique run identifier
SAMPLE1=sample1                   # Sample IDs for multi-individual MSMC analysis
SAMPLE2=sample2
SAMPLE3=sample3
SAMPLE4=sample4
CPIPSCRIPTS=/path/to/scripts      # Directory containing MSMC run scripts
POP1=pop1                         # Population identifier (group 1)
POP2=pop2                         # Population identifier (group 2)
Ne_POP1=10000                     # Effective population size for population 1
Ne_POP2=10000                     # Effective population size for population 2
N_boot=100                        # Number of bootstrap replicates

# Create output directory for the MSMC run
mkdir -p ${RUN_NAME}

# Generate multihetsep files for each chromosome (example: chromosomes 1, 2, and 3)
for CHR in 1 2 3; do
    # Accessible bases BED file for this chromosome
    REGION_TO_INCLUDE_perCHR_BED=/path/to/accessible_bases.chr${CHR}.bed

    # Mask files for each sample (assumed to be in the "msmc-prep" directory)
    MASKFILE_SAMPLE1=${SAMPLE1}.mask.chr${CHR}.bed.gz
    MASKFILE_SAMPLE2=${SAMPLE2}.mask.chr${CHR}.bed.gz
    MASKFILE_SAMPLE3=${SAMPLE3}.mask.chr${CHR}.bed.gz
    MASKFILE_SAMPLE4=${SAMPLE4}.mask.chr${CHR}.bed.gz

    # Final VCF files from MSMC prep for each sample
    MSMC_PREP_FINALVCF1=${SAMPLE1}.chr${CHR}.final.vcf.gz
    MSMC_PREP_FINALVCF2=${SAMPLE2}.chr${CHR}.final.vcf.gz
    MSMC_PREP_FINALVCF3=${SAMPLE3}.chr${CHR}.final.vcf.gz
    MSMC_PREP_FINALVCF4=${SAMPLE4}.chr${CHR}.final.vcf.gz

    /path/to/msmc-tools/generate_multihetsep.py \
        --chr NC_05186${CHR}.1 \
        --mask ${REGION_TO_INCLUDE_perCHR_BED} \
        --mask msmc-prep/${MASKFILE_SAMPLE1} \
        --mask msmc-prep/${MASKFILE_SAMPLE2} \
        --mask msmc-prep/${MASKFILE_SAMPLE3} \
        --mask msmc-prep/${MASKFILE_SAMPLE4} \
        msmc-prep/${MSMC_PREP_FINALVCF1} \
        msmc-prep/${MSMC_PREP_FINALVCF2} \
        msmc-prep/${MSMC_PREP_FINALVCF3} \
        msmc-prep/${MSMC_PREP_FINALVCF4} \
        > ${RUN_NAME}/${RUN_NAME}.chr${CHR}.multihetsep.txt &
done

wait

### Main MSMC2 run:
# Within-population coalescent analyses for two populations (2 individuals each)
msmc2 -t 4 -s -I 0-1,0-2,0-3,1-2,1-3,2-3 -o ${RUN_NAME}/${POP1} ${RUN_NAME}/*multihetsep*.txt &
msmc2 -t 4 -s -I 4-5,4-6,4-7,5-6,5-7,6-7 -o ${RUN_NAME}/${POP2} ${RUN_NAME}/*multihetsep*.txt &
# Cross-population coalescent analysis
msmc2 -t 8 -s -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 \
    -o ${RUN_NAME}/cross ${RUN_NAME}/*multihetsep*.txt &

# Combine and analyze MSMC2 results
/path/to/msmc-tools/combineCrossCoal.py \
    ${RUN_NAME}/cross.final.txt ${RUN_NAME}/${POP1}.final.txt ${RUN_NAME}/${POP2}.final.txt \
    > ${RUN_NAME}/${RUN_NAME}.combined.final.txt

## Bootstrapping: Create bootstrapped multihetsep files, which one can run msmc2 on
# /path/to/msmc-tools/multihetsep_bootstrap.py \
#     -n ${N_boot} \
#     -s 10000000 \
#     --chunks_per_chromosome 20 \
#     --nr_chromosomes 3 \
#     ${RUN_NAME}/boot \
#     ${RUN_NAME}/${RUN_NAME}.chr1.multihetsep.txt \
#     ${RUN_NAME}/${RUN_NAME}.chr2.multihetsep.txt \
#     ${RUN_NAME}/${RUN_NAME}.chr3.multihetsep.txt