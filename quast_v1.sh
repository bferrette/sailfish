#!/bin/bash
#
# Quast - Quality Assessment Tool for Genome Assemblies
# https://github.com/ablab/quast
#
# create conda environment
# conda creat -n QUASTenv -c bioconda -c conda-forge quast
#
# activate conda environment
# conda activate QUASTenv

# Run quast
quast.py \
        -o ./quast_results \
        --min-contig 1 \
        --thread 16 \
        --split-scaffolds \
        --large \
        --contig-thresholds 1000,5000,10000,50000,100000,500000,1000000,5000000,10000000 \
        *.fasta

