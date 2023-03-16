#!/bin/bash
#
# Merqury: k-mer based assembly evaluation
# https://github.com/marbl/merqury
#
# create conda environemt
conda create -n MERQURYenv -c conda-forge -c bioconda merqury openjdk=11

# Activate conda environment
conda activate MERQURYenv

# Best K-mer number
sh ${HOME}/anaconda3/envs/MERQURYenv/share/merqury/best_k.sh 619036510
genome: 619036510
tolerable collision rate: 0.001
19.5849

k=19.5849
for sample in ${HOME}/sailfish_cl/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_1.fq.gz ${HOME}/sailfish_cl/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_2.fq.gz
do
# 1. Build meryl dbs
    meryl k=$k count output $sample.meryl $sample threads=24 memory=48g
done

# 2. Merge
meryl union-sum output sailfish.meryl ./SFA29_FDSW202536543-1r_HKFGWDSXY_L4_*.meryl

# Run merqury
merqury.sh ./sailfish.meryl ${HOME}/sailfish_cl/assembly_cl/sailfish_CL_gap_closed_2.fasta sailfish.merqury.out
