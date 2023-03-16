#!/bin/bash
#
# Usage:
#   trim_reads_v2.sh <indir> <outdir_fastqs> <outdir_report> <njobs>
#
# Description:
#   Fastp: An ultra-fast all-in-one FASTQ preprocessor (QC/adapters/trimming/filtering/splitting/merging...)
#
# Requirements:
#   fastp - https://github.com/OpenGene/fastp
#   parallel - https://www.gnu.org/software/parallel/
#   multiqc - https://github.com/ewels/MultiQC
#
# Important:
#   Each job uses 8 CPU threads.

# trim reads and generate quality reports of before and after trimming in parallel
parallel \
  -j "$4" \
  --link \
  --plus \
  fastp \
    --in1 '{1}' \
    --in2 '{2}' \
    --out1 "$2"/'{1/}' \
    --out2 "$2"/'{2/}' \
    --detect_adapter_for_pe \
    --cut_tail \
    --cut_tail_window_size 4 \
    --cut_tail_mean_quality 15 \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --length_required 36 \
    --low_complexity_filter \
    --correction \
    --overrepresentation_analysis \
    --json "$3"/'{1/..}-2.fastp.json' \
    --html "$3"/'{1/..}-2.fastp.html' \
    --report_title '{1/..}-2' \
    --thread 8 \
    --compression 9 \
  ::: "$1"/*_1.fq.gz \
  ::: "$1"/*_2.fq.gz

# Compile all reports with MultiQC
multiqc . --interactive
