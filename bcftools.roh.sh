#!/bin/bash
#
# Detecting runs of homozygosity (RoH) - https://samtools.github.io/bcftools/howtos/roh-calling.html
#
# Description:
# BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data
# https://doi.org/10.1093/bioinformatics/btw044
#
# Usage:
# sh bcftools.roh.sh <indir> <outdir> <ncores> <njobs>
# 
# bcftools RoH cluster atlantic
parallel -j $4 bcftools roh --threads $3 --skip-indels -r {} --estimate-AF "PL,ATL.cluster" -S ATL.cluster $1/snps.filtered.bcf > ATL.roh.txt ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32
for f in *.txt; do sample=$(basename $f .roh.txt); grep '^RG' $f | tail -n +2 > ATL.roh.tsv; done

# bcftools RoH cluster indo-pacific
parallel -j $4 bcftools roh --threads $3 --skip-indels -r {} --estimate-AF "PL,IDWP.cluster" -S IDWP.cluster $1/snps.filtered.bcf > IDWP.roh.txt ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32
for f in *.txt; do sample=$(basename $f .roh.txt); grep '^RG' $f | tail -n +2 > IDWP.roh.tsv; done

# compress all intermediary files
pigz --best *.txt
