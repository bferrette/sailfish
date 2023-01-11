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

# bcftools RoH cluster atlantic per sample
for sample in SFA01 SFA02 SFA03 SFA04 SFA05 SFA06 SFA07 SFA08 SFA09 SFA10 SFA11 SFA12 SFA13 SFA14 SFA15 SFA16 SFA18 SFA17 SFA19 SFA20 SFA21 SFA22 SFA23 SFA24 SFA25 SFA26 SFA27 SFA28 SFA29 SFA30 SFA31 SFA32; do
	parallel -j $4 bcftools roh --threads $3 --skip-indels -r {} --estimate-AF "PL,ATL.cluster" -s ${sample} $1/snps.filtered.bcf > ${sample}.roh.txt ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep '^RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# bcftools RoH cluster indo-pacific per sample
for sample in SFA33 SFA34 SFA35 SFA36 SFA37 SFA38 SFA39 SFA40 SFA41 SFA42 SFA43 SFA44 SFA45 SFA46 SFA47 SFA48 SFA49 SFA50 SFA51 SFA52 SFA53 SFA54 SFA55 SFA56 SFA57 SFA58 SFA59 SFA60 SHS-SA01; do
	parallel -j $4 bcftools roh --threads $3 --skip-indels -r {} --estimate-AF "PL,IDWP.cluster" -s ${sample} $1/snps.filtered.bcf > ${sample}.roh.txt ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# compress all intermediary files
pigz --best *.txt

