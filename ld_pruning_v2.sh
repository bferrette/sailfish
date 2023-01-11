#!/bin/bash
#
# Usage:
#   ld_pruning_v2.sh <indir> <outdir> <ncpus>
#
# Description:
#   Calculate pairwise linkage disequilibrium with ngsLD, plot LD decay curve
#   with 'fit_LDdecay.R', prune linked sites with 'prune_graph.pl', and extract
#   unlinked sites from ANGSD's genotype likelihoods file.
#
# Requirements:
#   parallel
#   ngsld (incl. fit_LDdecay.R, prune_graph.pl)
#   r

# path to ngsLD
ngsld=${HOME}/software/ngsLD

# number of samples
n_ind=$(cat $1/bamlist | wc -l)
# find MAF file
maf=$(find $1 -name '*.mafs.gz')
# create file with SNP positions
zcat ${maf} | cut -f 1,2 | tail -n +2 > $2/snps.pos
# count the number of SNPs
n_sites=$(cat $2/snps.pos | wc -l)
# find genotype likelihoods (GL) file
gl=$(find $1 -name '*.beagle.gz')

# calculate linkage disequilibrium (LD) for all SNP pairs up to 1 Mb apart
${ngsld}/ngsLD \
  --geno ${gl} \
  --probs \
  --n_ind ${n_ind} \
  --n_sites ${n_sites} \
  --pos $2/snps.pos \
  --max_kb_dist 200 \
  --min_maf 0.05 \
  --out $2/snps.ld \
  --n_threads $3 \
  &> $2/ngsld.log

# sample ngsLD output for fitting LD decay curve
awk 'rand()<0.0001' $2/snps.ld > $2/snps.sampled.ld

# fit LD decay curve
ls $2/snps.sampled.ld \
  | Rscript \
    --vanilla \
    --slave ${ngsld}/scripts/fit_LDdecay.R \
    --ld r2 \
    --n_ind ${n_ind} \
    --max_kb_dist 100 \
    --fit_boot 100 \
    --fit_bin_size 250 \
    --fit_level 100 \
    --plot_data \
    --plot_scale 3 \
    --out $2/ld_decay.pdf \
    &> $2/fit_LDdecay.log

# START HERE
# split ngsLD output to separate files by scaffold
parallel -j $3 "grep -P {}: $2/snps.ld > $2/snps.ld.{}" ::: chr{1..14}

# prune linked sites by scaffolds in parallel
parallel -j $3 \
  "${ngsld}/scripts/prune_graph.pl --in_file {} --max_kb_dist 25 --min_weight 0.1 --out {}.pruned &> {}.prune_graph.log" \
  ::: $2/snps.ld.*

# concatenate and sort unlinked sites into a single file
cat $(ls $2/snps.ld.chr*.pruned) | sort -V | sed 's/:/_/' > $2/snps.ld.pruned
# concatenate 'prune_graph.pl' log files
cat $(ls -v $2/snps.ld.chr*.prune_graph.log) > $2/prune_graph.log
# remove intermediate files
rm $2/snps.ld.chr*

# set output file name
pruned_gl=$2/snps.ld_pruned.beagle
# extract SNPs matching the '*.ld.pruned' file from the GL file
zcat ${gl} | grep -F -f $2/snps.ld.pruned > ${pruned_gl}.tmp
# get linked SNPs inadvertently extracted due to partial string matching
awk '{ print $1 }' ${pruned_gl}.tmp \
  | diff $2/snps.ld.pruned - \
  | grep -Eo 'chr.+' \
  > $2/sites2remove
# extract the header of the original GL file
zcat ${gl} | head -1 > ${pruned_gl}
# extract only the unlinked SNPs from the GL (.tmp) file
grep -v -F -f $2/sites2remove ${pruned_gl}.tmp >> ${pruned_gl}
# compress the pruned GL file
gzip ${pruned_gl}
# remove intermediate files
rm ${pruned_gl}.tmp
