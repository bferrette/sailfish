#!/bin/bash
#
# Usage:
#   genotype_calling_v2.sh <ref.fa> <indir_bams> <site_depth.stats> <outdir> <nthreads> <njobs>
#
# Description:
#   Call genotypes with bcftools.
#
# Requirements:
#   bcftools

# create a list of BAM files sorted by population
ls -1 -v $2/*.clean.bam >> $4/bamlist

# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' $3)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $3)
# if minimum depth threshold is negative make it 1
if [ ${min_dp} -lt 0 ]; then min_dp=1; fi

# joint genotype calling
parallel -j $6 --tmpdir ${HOME}/tmp \
"bcftools mpileup --threads $5 -b $4/bamlist -C 50 -f $1 -r {} -q 30 -Q 30 --rf 2 --ff 4,256,512,1024 -Ov | bcftools call --threads $5 -mv -a GQ,GP -Ov -o $4/{}.vcf" ::: scaff{1..16896}

# Concatenate all vcf of chromosomes:
bcftools concat --threads $5 $4/scaff*.vcf -Ob -o $4/snps.bcf && rm $4/scaff*.vcf

# filter multisample BCF
#bcftools filter --threads $5 -e "DP<${min_dp} || DP>${max_dp} || MQ<30 || QUAL<30 || F_MISSING>=0.05 || MAF[0]<0.05" -g 3 -s 'LowQual' -Ou $4/snps.bcf \
#  | bcftools view --threads $5 -m 2 -M 2 -v snps -c 1:minor -i 'FILTER="PASS"' -Ob -o $4/snps.filtered.bcf -

bcftools filter --threads $5 -e "DP<${min_dp} || DP>${max_dp} || MQ<30 || QUAL<30" -g 3 -s 'LowQual' -Ou $4/snps.bcf \
  | bcftools view --threads $5 -m 2 -M 2 -v snps -c 1:minor -i 'FILTER="PASS"' -Ob -o $4/snps.filtered.bcf -

# index BCF
bcftools index $4/snps.filtered.bcf
