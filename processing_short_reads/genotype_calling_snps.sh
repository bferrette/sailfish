#!/bin/bash
#
# Usage:
#   genotype_calling.sh <ref.fa> <indir_bams> <site_depth.stats> <outdir> <ncores> <njobs>
#
# Description:
#   Call genotypes with bcftools.
#
# Requirements:
#   bcftools - https://github.com/samtools/bcftools
#   gnu parallel - https://www.gnu.org/software/parallel/

# create a list of BAM files sorted by population
ls -1 -v $2/*.clean.bam >> $4/bamlist

# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' $3)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $3)
# if minimum depth threshold is negative make it 1
if [ ${min_dp} -lt 0 ]; then min_dp=1; fi

# parallel genotype calling and sort all variable sites per chromosome of clean BAM files
parallel -j $6 --tmpdir /home/bferrette/tmp \
"bcftools mpileup --threads $5 -b $4/bamlist -C 50 -f $1 -r {} -q 30 -Q 30 --rf 2 --ff 4,256,512,1024 -Ou | bcftools call --threads $5 -mv -a GQ,GP -Ou | bcftools sort -Ob -o $4/{}.sorted.bcf" ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32

# All sites
#bcftools mpileup --threads $5 -b $4/bamlist -C 50 -f $1 -q 30 -Q 30 --rf 2 --ff 4,256,512,1024 -Ou | bcftools call --threads $5 -mv -a GQ,GP | bcftools sort -Ob -o $4/snps.sorted.bcf -
# ANGSD sites
#bcftools mpileup --threads $5 -b $4/bamlist -C 50 -f $1 -T sites.txt -q 30 -Q 30 --rf 2 --ff 4,256,512,1024 -Ou | bcftools call --threads $5 -mv -f GQ,GP -Ob -o $4/snps.sorted.bcf -

# index all sorted VCFs
parallel -j $6 "bcftools index $4/{}" ::: *.sorted.bcf

# Concatenate all sorted VCFs of all chromosomes:
bcftools concat --threads $5 $4/*.sorted.bcf -Ov -o $4/snps.sorted.bcf

# filter multisample BCF
bcftools filter --threads $5 -e "DP<=${min_dp} || DP>=${max_dp} || MQ<=30 || QUAL<=30" -g 3 -s 'LowQual' -Ou $4/snps.sorted.bcf \
  | bcftools view --threads $5 -m 2 -M 2 -v snps -c 1:minor -i 'FILTER="PASS"' -Ob -o $4/snps.filtered.sorted.bcf -

# filter multisample VCF by genotype quality of SNPs and include only biallelic SNPs ***MORE STRICT FILTERING***
#bcftools filter --threads $5 -e "DP<${min_dp} || DP>${max_dp} || GQ<20" -g 3 -s 'LowQual' -Ou $4/snps.sorted.bcf \
#  | bcftools view --threads $5 -m 2 -M 2 -v snps -c 1:minor -i 'FILTER="PASS"' -Ob -o $4/snps.filtered.sorted.bcf -

# index BCF
bcftools index $4/snps.filtered.sorted.bcf

# extract all filtered SNPs per chromosome from the filtered multisample VCF
parallel -j $6 --tmpdir /home/bferrette/tmp \
"bcftools convert --threads $5 -r {} $4/snps.filtered.sorted.bcf -Oz -o $4/{}.filtered.sorted.vcf.gz" ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32
# index all sorted VCFs
parallel -j $6 "bcftools index $4/{}" ::: *.vcf.gz

