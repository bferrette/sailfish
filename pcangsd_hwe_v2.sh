#!/bin/bash
#
# Usage:
#   pcangsd_hwe.sh <indir> <outdir> <ncpus>
#
# Description:
#   Calculate a covariance matrix for SNPs filtering out those that strongly
#   deviate from the Hardy-Weinberg equilibrium using PCAngsd.
#
# Requirements:
#   pcangsd
#   python3

# path to PCAngsd
pcangsd=${HOME}/software/pcangsd/pcangsd.py

# find genotype likelihoods (GL) file
gl=$(find $1 -name '*.beagle.gz')
# set output name and directory
out=$2/$(basename ${gl%.beagle.gz})

# calculate covariance matrix, estimate per-site inbreeding coefficients, and
# perform likehood ratio tests for HWE
python3 ${pcangsd} \
  -beagle ${gl} \
  -minMaf 0.05 \
  -o ${out} \
  -inbreedSites \
  -threads $3 \
  &> ${out}.pcangsd.log

# add delimiter between logs to 'pcangsd.log'
echo -e "\n---\n" >> ${out}.pcangsd.log

# calculate covariance matrix with HWE filter
python3 ${pcangsd} \
  -beagle ${gl} \
  -minMaf 0.05 \
  -o ${out}.hwe_filter \
  -hwe ${out}.lrt.sites.npy \
  -hwe_tole 1e-6 \
  -threads $3 \
  &>> ${out}.pcangsd.log
