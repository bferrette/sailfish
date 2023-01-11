#!/bin/bash
#
# Usage:
#   realSFS_per_pop.sh <ref.fa> <anc.fa> <indir_bam> <pop.site_depth.stats> <outdir>
#
# Description:
#   Estimate per sample site frequency spectrum per population with ANGSD.
#
# Requirements:
#   angsd
#   parallel
#
# Important:
#   Each job uses up to 4 CPU threads.

# create a list of BAM files sorted by population
ls -1 -v $3/*.clean.bam >> $5/bamlist

# count the number of individuals
n_ind=$(cat $5/bamlist | wc -l)
# set minimum number of individuals per site
min_ind=$(python -c "print(f'{round(${n_ind} * 0.9)}')")
# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' $4)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $4)
# if minimum depth threshold is negative make it 1
if [[ ${min_dp} -lt 0 ]]; then min_dp=1; fi
# set filters
filters="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -minQ 30 -setMaxDepth ${max_dp} -setMinDepth ${min_dp} -skipTriallelic -sb_pval 1e-6 -hetbias_pval 1e-6"
# set tasks
todo="-doSnpStat 1 -doHWE 1 -doCounts 1 -doMajorMinor 1 -doMaf 2 -doPost 1 -doGeno 8 -doSaf 1"
# estimate the unfolded sample allele frequency (SAF) likelihoods for each site 
angsd -b $5/bamlist -ref $1 -anc $2 ${filters} ${todo} -GL 1 -P 4 -out $5/pop &> $5/pop.angsd.log
# estimate the folded site frequency spectrum (SFS)
#realSFS $5/pop.saf.idx -fold 1 > $5/pop.sfs 2> $5/pop.realSFS.log
# estimate the unfolded site frequency spectrum (SFS)
realSFS $5/pop.saf.idx > $5/pop.sfs 2> $5/pop.realSFS.log
