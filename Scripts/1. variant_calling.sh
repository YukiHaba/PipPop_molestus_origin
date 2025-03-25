#!/bin/bash

##############################
#### Raw Variant Calling ####
##############################
# Variables used in this section:
#   BAM_LIST: File containing a list of BAM files.
#   REFERENCE_GENOME: Path to the reference genome FASTA file. CpipJ5 assembly (NCBI RefSeq assembly: GCF_015732765.1)
#   THREADS: Number of threads for parallel processing.
#   REGION: Genomic region to process. 20Mb chunks (e.g., "NC_051861.1:1-20000000").
#   HEADER: Prefix for naming output files.

bcftools mpileup \
    -b ${BAM_LIST} \                                                
    -f GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna \    
    --threads ${THREADS} \                                          
    -r ${REGION} \                                                  
    -Ou -a FORMAT/DP,FORMAT/AD,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR | \
bcftools call \
    --threads ${THREADS} \                          
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

# Variables used in this section:


# run this per region (e.g., "NC_051861.1:1-20000000")
bcftools filter -g 5 $ALLSITES_VARIANTS -Ou | # mask 5bps around indels
bcftools view -m2 -M2 -v snps -Ou | # get only biallelic snps
bcftools view -Ob -T ^$REPEATBED \  # GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_rm.out.bed 
-i 'F_MISSING < 0.25 & MQ > 40 & QUAL > 30' \
> /tmp/$REGION.good_biallelic_snps.rm.bcf

bcftools index /tmp/$REGION.good_biallelic_snps.rm.bcf

# Concatenate the REGION filtered VCF files into a single file.
bcftools concat \
    --allow-overlaps \
    --file-list ${VARIANT_FILE_LIST} \   # Text file listing all regional VCF files.
    --threads 20 \                       # Use 20 threads for the operation.
    -Oz -o /tmp/good_biallelic_snps.rm.bcf_combined.vcf.gz  # Output concatenated VCF file.
    
# Index the concatenated VCF.
bcftools index /tmp/good_biallelic_snps.rm.bcf_combined.vcf.gz 

# Extract variants that fall in accessible genomic regions.
# Variables:
#   ACCESSIBLE_BED: BED file with accessible regions for the chromosome. data/PipPop.accessible_sites.bed
#   INPUT_BCF: Name of the input BCF file containing filtered variants.
ACCESSIBLE_BED="/path/to/accessible_regions/PipPop.accessible_sites.bed"  # Accessible regions BED file.
INPUT_VCF="good_biallelic_snps.rm.bcf_combined.vcf.gz"     # Adjust input file naming as needed.

bcftools view \
    ${INPUT_VCF} \                       # Input BCF file with filtered variants.
    -R ${ACCESSIBLE_BED} \               # Restrict to accessible regions.
    -Oz -o good_biallelic_snps.rm.bcf_combined.accessible.vcf.gz &  # Output VCF with accessible variants.


##############################
# (for some downstream analysis) get unlinked SNPs
##############################
VCFHEADER=good_biallelic_snps.rm.bcf_combined.accessible

## 1. vcf to map/ped
plink --vcf $VCFHEADER.vcf.gz \
--allow-extra-chr \
--recode --out /tmp/$VCFHEADER \
--set-missing-var-ids @:#[Cpip] \
--threads 20

## 2. to get a list of SNPs R2 >= 0.2
## take chunks of 200 snps, sliding 20 snps, remove all pairs being r^2 > 0.2
plink --file /tmp/$VCFHEADER \
--allow-extra-chr \
--indep-pairwise 200 20 0.2 \
--out /tmp/$VCFHEADER \
--threads 20

## 3. perform the extraction of the SNPs and convert to vcf
plink --file /tmp/$VCFHEADER \
--allow-extra-chr \
--extract /tmp/$VCFHEADER.prune.in \
--out /tmp/$VCFHEADER.LDpruned \
--recode vcf \
--threads 20

## 4. bgzip and index
bgzip /tmp/$VCFHEADER.LDpruned.vcf
# rename samples in the file
bcftools reheader -s /tmp/$VCFHEADER.samplenames /tmp/$VCFHEADER.LDpruned.vcf -o /tmp/$VCFHEADER.LDpruned.renamed.vcf
rm /tmp/$VCFHEADER.LDpruned.vcf
mv /tmp/$VCFHEADER.LDpruned.renamed.vcf /tmp/$VCFHEADER.LDpruned.vcf
# index
bcftools index /tmp/$VCFHEADER.LDpruned.vcf.gz
