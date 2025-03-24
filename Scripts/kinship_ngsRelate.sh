#!/bin/bash

# Variables
THREADS=20                                     
N_IND=840                                                        

# Run NGSrelate to compute kinship coefficients:
ngsRelate \
    -h PipPop_all.840ind.vcf.gz \             # input gzipped VCF file
    -z ind_name_file.txt \                    # file with individual names
    -p ${THREADS} \                           # Set the number of processing threads
    -n ${N_IND} \                             # Number of individuals in the study
    -O /tmp/ngsRelate.${OUT_HEADER}           # Define the output file prefix (results saved in /tmp)