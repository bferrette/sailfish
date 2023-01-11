#!/bin/bash
#
# Usage:
#   psmc_v2.sh <ref.fa> <indir> <outdir> <dp_persample_dir> <ncpus> <njobs>
#
# Description:
#   This software package infers population size history from a diploid sequence
#   using the Pairwise Sequentially Markovian Coalescent (PSMC) model.
#   (https://github.com/lh3/psmc)
#
# Requirements:
#   psmc
#   parallel
#
samples=(SFA01 SFA02 SFA03 SFA04 SFA05 SFA06 SFA07 SFA08 SFA09 SFA10 SFA11 SFA12 SFA13 SFA14 SFA15 SFA16 SFA17 SFA18 SFA19 SFA20 SFA21 SFA22 SFA23 SFA24 SFA25 SFA26 SFA27 SFA28 SFA29 SFA30 SFA31 SFA32 SFA33 SFA34 SFA35 SFA36 SFA37 SFA38 SFA39 SFA40 SFA41 SFA42 SFA43 SFA44 SFA45 SFA46 SFA47 SFA48 SFA49 SFA50 SFA51 SFA52 SFA53 SFA54 SFA55 SFA56 SFA57 SFA58 SFA59 SFA60 SFA61.20x)
#ref=$HOME/sailfish_cl/assembly_cl/sailfish_CL_gap_closed_2.fasta
#dp_stats=$HOME/sailfish_cl/processing/depth_stats/per_sample
tmp=$HOME/tmp
# run PSMC
for sample in ${samples[@]}; do
  median_dp=$(grep -Po 'Median: \K\d+' $4/${sample}*.depth.stats)
  max_dp=$(python3 -c "print(f'{int($median_dp) * 2}')")
  min_dp=$(python3 -c "print(f'{int($median_dp)/3}')")
  echo "bcftools mpileup --threads $5 -C 50 -f $1 -q 30 -Q 30 --rf 2 --ff 4,256,512,1024 -Ou $2/${sample}.clean.bam | \
        bcftools call --threads $5 -c - | \
        vcfutils.pl vcf2fq -d ${min_dp} -D ${max_dp} -Q 30 | pigz > ${sample}.consensus.fq.gz" >> consensus_call.jobs
  echo "fq2psmcfa -q 30 ${sample}.consensus.fq.gz > ${sample}.psmcfa" >> fq2psmcfa.jobs
  echo "splitfa ${sample}.psmcfa 10000 > ${sample}.split.psmcfa" >> splitfa.jobs
  echo "psmc -N 25 -t 15 -r 5 -p '4+25*2+4+6' -o ${sample}.psmc ${sample}.psmcfa" >> psmc.jobs
  for round in {1..100}; do
    echo "psmc -N 25 -t 15 -r 5 -b -p '4+25*2+4+6' -o ${sample}.round-$(printf '%03d' $round).psmc ${sample}.split.psmcfa" >> psmc_boot.jobs
  done
done

cat consensus_call.jobs | parallel -j $6 --tmpdir ${tmp} && \
cat fq2psmcfa.jobs | parallel -j $6 --tmpdir ${tmp} && \
cat splitfa.jobs | parallel -j $6 --tmpdir ${tmp} && \
cat psmc.jobs | parallel -j $6 --tmpdir ${tmp} && \
cat psmc_boot.jobs | parallel -j $6 --tmpdir ${tmp}
