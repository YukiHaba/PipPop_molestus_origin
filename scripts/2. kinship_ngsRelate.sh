#!/bin/bash

# Variables
THREADS=20                                     
N_IND=840                                                        

# Run NGSrelate to compute kinship coefficients:
ngsRelate \
    -h good_biallelic_snps.rm.bcf_combined.accessible.vcf.gz \
    -z ind_name_file.txt \                    # file with individual names
    -p ${THREADS} \                           # Set the number of processing threads
    -n ${N_IND} \                             # Number of individuals in the study
    -O /tmp/ngsRelate.stats                   # Define the output file prefix (results saved in /tmp)