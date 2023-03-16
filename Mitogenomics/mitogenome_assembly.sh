#!/bin/bash
# Mitogenomes assembly
#
# Usage: mitogenome_assembly.sh <indir> <outdir> <ncores>
#
# GetOrganelle - https://github.com/Kinggerm/GetOrganelle
# This toolkit assemblies organelle genome from genomic skimming data.

# install GetOrganelle
# conda create -n GETORGANELLEenv -c bioconda getorganelle blast=2.13 perl=5.26 bowtie2=2.4.5

# activate conda environment
conda activate GETORGANELLEenv

# After installation of GetOrganelle v1.7.0+, download and initialize the database of your organelle genome
get_organelle_config.py --add animal_mt

# iterate over read 1 files
for r1 in "$1"/*_1.fq.gz; do
# read 2
r2="${r1/_1./_2.}"
# basename
base="$(basename "${r1/_1.fq.gz}")"
# run to assembly animal mitogenomes
get_organelle_from_reads.py -1 $r1 -2 $r2 -R 10 -k 21,45,65,85,105 -F animal_mt -o $2/${base} -t $3
done

