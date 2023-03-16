#!/bin/bash
#
# Usage:
#   clean_bams_v2.sh <regions.bed> <indir> <outdir> <nthreads ><njobs>
#
# Description:
#   Remove unmapped (4), secondary (256), QC failed (512), duplicate (1024), and
#   supplementary (2048) reads from indel-realigned BAMs, and keep only reads
#   mapped in a proper pair (2) to regions in a BED file (non-repetitive regions
#   in autosomes) using SAMtools.
#
# Requirements:
#   samtools - https://github.com/samtools/samtools
#   gnu parallel - https://www.gnu.org/software/parallel/
#
# clean indel-realigned BAMs
parallel -j "$5" --plus \
  "samtools view -@ $4 -b -F 3844 -f 2 -L $1 -o $3/{/..}.clean.bam {}" \
  ::: "$2"/*.realigned.bam

# index clean BAMs
parallel -j "$5" samtools index -b {} ::: "$3"/*.clean.bam

