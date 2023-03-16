#!/bin/bash
#
# Usage:
#   snp_calling.sh <ref.fa> <indir_bams> <site_depth.global.stats> <outdir>
#
# Description:
#   Call and filter SNPs with ANGSD.
#
# Requirements:
#   angsd - https://github.com/ANGSD/angsd
#   python3 - https://www.python.org/
#
# create a list of BAM files sorted by population
ls -1 -v $2/*.clean.bam >> $4/bamlist

# count the number of individuals
n_ind=$(cat $4/bamlist | wc -l)
# set minimum number of individuals per site
min_ind=$(python3 -c "print(f'{round(${n_ind} * 0.9)}')")
# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' $3)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $3)
# if minimum depth threshold is negative make it 1
if [[ ${min_dp} -lt 0 ]]; then min_dp=1; fi

# SNP calling with filters
/opt/software/angsd/angsd -b $4/bamlist -ref $1 -GL 1 -P 8 -out $4/sailfish \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 -baq 1 \
  -minMapQ 30 -minQ 30 -minInd ${min_ind} \
  -doSnpStat 1 -doHWE 1 -sb_pval 1e-6 -hwe_pval 1e-6 -hetbias_pval 1e-6 \
  -doCounts 1 -setMinDepth ${min_dp} -setMaxDepth ${max_dp} \
  -doMajorMinor 1 -skipTriallelic 1 \
  -doMaf 2 -doPost 1 -minMaf 0 -SNP_pval 1e-6 \
  -doGeno 8 -geno_minDepth 3 -doGlf 2 \
  &> $4/angsd.log

