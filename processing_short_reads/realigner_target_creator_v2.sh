#!/bin/bash
#
# Usage:
#   realigner_target_creator_v2.sh <ref.fa> <indir> <outdir> <ndata_threads> <njobs>
#
# Description:
#   Identifies regions where alignments may potentially be improved using
#   GATK RealignerTargetCreator.
#
# Requirements:
#   gatk 3.8-1 - https://github.com/broadinstitute/gatk/releases
#   opejdk 1.8 - https://openjdk.org/projects/jdk8u/
#   gnu parallel - https://www.gnu.org/software/parallel/
#   picard - https://github.com/broadinstitute/picard
#   samtools - https://github.com/samtools/samtools
#
# Important:
#   Each job will use 'N data threads' + 4 CPUs ('-nt' + ParallelGCThreads).
#   Each data thread will use 16GB RAM.
#
# path to Picard
picard=${HOME}/software/picard/picard.jar

# path to GATK v3.8-1
gatk=${HOME}/software/GenomeAnalysisTK.jar

# index reference genome FASTA
samtools faidx $1

# create a sequence dictionary for the reference genome FASTA
#java -jar ${picard} CreateSequenceDictionary R=$1 O=${1/.fasta/.dict}

# get formatted string with input BAMs
bams=$(find $2 -name '*.dedup.bam' -printf '-I %p ')

# get chromosome IDs and spans formatted as 'chr:start-end'
cat $1.fai | awk '{ print $1":1-"$2}' \
  | while read chr; do # iterate over chromosomes
    # echo command to a file
    echo \
      "/usr/bin/java -Xmx16G -XX:ParallelGCThreads=4 -jar ${gatk} \
        -T RealignerTargetCreator \
        -nt $4 \
        -R $1 \
        -L ${chr} \
        ${bams} \
        -o $3/${chr%:*}.intervals \
        &> $3/${chr%:*}.rtc.log" \
      >> $3/rtc.jobs
  done

# run RealignerTargetCreator in parallel
cat $3/rtc.jobs | parallel -j $5

# search for warn messages
n_warns=$(grep -L 'Done. There were no warn messages.' $3/*.rtc.log | wc -l)

# if there are no warnings, concatenate intervals and remove intermediate files
if [[ ${n_warns} -eq 0 ]]; then
  ls -v $3/*.intervals > $3/chr.intervals.list &&
  cat $(cat $3/chr.intervals.list) > $3/realigner.intervals &&
  echo 'RealignmentTargetCreator completed successfully!' &&
  rm $3/*.{jobs,list}

# if there are warnings, print a list of files in which they occurred
else
  echo 'These files had warn messages:'
  grep -L 'Done. There were no warn messages.' $3/*.rtc.log
fi

