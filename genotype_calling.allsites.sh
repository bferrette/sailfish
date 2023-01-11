#!/bin/bash
#
# Usage:
#   genotype_calling.sh <ref.fa> <indir_bams> <site_depth.stats> <outdir> <ncores> <njobs>
#
# Description:
#   Call genotypes with bcftools.
#
# Requirements:
#   bcftools

# create a list of BAM files sorted by population
#ls -1 -v $2/*.clean.bam >> $4/bamlist

# set minimum and maximum depth thresholds
min_dp=$(grep -Po '\(MEDIAN \- 5 \* MAD\): \K\-*\d+' $3)
max_dp=$(grep -Po '\(MEDIAN \+ 5 \* MAD\): \K\d+' $3)
# if minimum depth threshold is negative make it 1
if [ ${min_dp} -lt 0 ]; then min_dp=1; fi

# joint genotype calling of all variable and invariable sites
parallel -j $6 --tmpdir /home/bferrette/sailfish_cl/population.analysis/sailfish/heterozygozity/allsites.allsamples/tmp \
"bcftools mpileup --threads $5 -b $4/bamlist -C 50 -f $1 -r {} -q 30 -Q 30 --rf 2 --ff 4,256,512,1024 -Ov | bcftools call --threads $5 -m -a AD,GQ,GP -Ov -o $4/{}.vcf" ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32

# Concatenate all vcf of chromosomes:
bcftools concat --threads $5 HiC_scaffold_*.vcf -Ov -o sailfish.vcf
# filter VCF (optional)
bcftools filter --threads $5 -e "DP<${min_dp} || DP>${max_dp} || MQ<30 || QUAL<30 || F_MISSING>=0.1" -g 3 -s 'LowQual' sailfish.vcf | bcftools view --threads $5 -i 'FILTER="PASS"' -Oz -o $4/sailfish.filtered.vcf.gz -

# compress all vcf files
parallel -j $6 --tmpdir /home/bferrette/sailfish_cl/population.analysis/sailfish/heterozygozity/allsites.allsamples/tmp \
"bcftools convert --threads 12 {}.vcf -Oz -o {}.vcf.gz" ::: HiC_scaffold_5 HiC_scaffold_6 HiC_scaffold_7 HiC_scaffold_9 HiC_scaffold_11 HiC_scaffold_12 HiC_scaffold_13 HiC_scaffold_16 HiC_scaffold_17 HiC_scaffold_18 HiC_scaffold_19 HiC_scaffold_20 HiC_scaffold_21 HiC_scaffold_22 HiC_scaffold_23 HiC_scaffold_24 HiC_scaffold_25 HiC_scaffold_26 HiC_scaffold_27 HiC_scaffold_28 HiC_scaffold_29 HiC_scaffold_30 HiC_scaffold_31 HiC_scaffold_32

# compress vcf file
bcftools convert --threads $5 sailfish.vcf -Oz -o sailfish.vcf.gz

# index all VCFs
parallel -j $6 "bcftools index $4/{}" ::: *.vcf.gz 
