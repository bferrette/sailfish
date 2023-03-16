#!/bin/bash
#
# Usage:
#   site_depth_global.sh <indir> <outdir> <ncpus> <njobs>
#
# Description:
#   Calculate the global site depth from cleaned BAMs using sambamba.
#
# Requirements:
#   gnu parallel - https://www.gnu.org/software/parallel/
#   sambamba - https://github.com/biod/sambamba

# set shell extglob
shopt -s extglob

# create an array of giraffe BAMs
bams=($1/*.clean.bam)

# unset shell extglob
shopt -u extglob

sambamba depth base -t $3 --combined --fix-mate-overlaps ${bams[@]} | awk 'NR>1 { print $3 }' > $2/site_depth.global

# randomly sample 0.01% of all sites for plotting the global site depth distribution
awk 'rand()<0.0005' $2/site_depth.global > $2/site_depth.global.sampled

# calculate summary statistics for the global site depth
python3 site_depth_stats.py site_depth.global > site_depth.global.stats

# calculte per sample depth stats
for bam in *.clean.bam; do
  sample=$(basename ${bam%.clean.bam})
  sambamba depth base -t $3 --fix-mate-overlaps ${bam} | awk 'NR>1 { print $3 }' > ${sample}.depth
  python3 site_depth_stats.py ${sample}.depth > ${sample}.depth.stats
done
