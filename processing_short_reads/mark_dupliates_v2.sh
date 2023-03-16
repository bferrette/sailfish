#!/bin/bash
#
# Usage:
#   mark_duplicates_v2.sh <indir> <outdir_bam> <outdir_log> <njobs>
#
# Description:
#   Mark PCR/optical duplicate reads with Picard's MarkDuplicates.
#
# Requirements:
#   gnu parallel - https://www.gnu.org/software/parallel/
#   picard - https://github.com/broadinstitute/picard
#   samtools - https://github.com/samtools/samtools
#   opejdk 1.8 - https://openjdk.org/projects/jdk8u/
#
# Important:
#   Each Picard job uses 4 CPU threads.

# path to Picard
picard=${HOME}/software/picard/build/libs/picard.jar

# mark PCR/optical duplicate reads for patterned flowcells
parallel \
  -j "$4" \
  --plus \
  "/urs/bin/java -XX:ParallelGCThreads=4 -Xmx16G -jar ${picard} \
    MarkDuplicates \
      I={} \
      O=$2/{/..}.dedup.bam \
      M=$3/{/..}.dedup.metrics.txt \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      &> $3/{/..}.markduplicates.log" \
  ::: "$1"/*.sorted.bam

# index new BAMs
parallel -j "$4" samtools index -b {} ::: "$2"/*.dedup.bam

