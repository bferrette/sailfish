#!/bin/bash
#
# Usage:
#   check_mapping_quality_v2.sh <indir> <outdir_flagstats> <outdir_qualimap> <ncpus>
#
# Description:
#   Calculate mapping statistics with SAMtools and generate plots of various
#   mapping quality measurements with QualiMap.
#
# Requirements:
#   parallel - https://www.gnu.org/software/parallel/
#   qualimap - http://qualimap.conesalab.org/
#   samtools - https://github.com/samtools/samtools
#
# calculate mapping statistics
parallel -j $4 "samtools flagstat {} > $2/{/.}.flagstats" ::: $1/*.dedup.bam

# iterate over elements in the array
for bam in $1/*.dedup.bam; do
  # output directory
  outdir=$3/$(basename ${bam/bam/qualimap})
  # generate mapping quality report
  qualimap bamqc \
    -bam ${bam} \
    --collect-overlap-pairs \
    --skip-duplicated \
    -nt $4 \
    -outdir ${outdir} \
    --java-mem-size=24G \
    &> ${outdir}.log
done
