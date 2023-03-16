#!/bin/bash
#
# BUSCOv5 - Benchmarking sets of Universal Single-Copy Orthologs
# https://gitlab.com/ezlab/busco/-/tree/master
#
# create conda environemt
# conda create -n BUSCOenv -c conda-forge -c bioconda busco=5.4.6
#
# activate conda environment
# conda activate BUSCOenv
#
# run BUSCO in genome mode
busco --cpu 32 --mode genome --long --in Istiophorus_platypterus.fasta --out busco --lineage vertebrata_odb10 --augustus --augustus_species zebrafish
