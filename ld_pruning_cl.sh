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
#   ngsLD (incl. fit_LDdecay.R, prune_graph.pl) - https://github.com/fgvieira/ngsLD
#   R

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
  --min_maf 0 \
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
parallel --tmpdir /home/bferrette/sailfish_cl/population.analysis/sailfish/ld_pruning2/tmp \
-j $3 "grep -P {}: $2/snps.ld > $2/snps.ld.{}" ::: $(cat /home/bferrette/sailfish_cl/assembly_cl/repetitive_regions/repeatmasker_step2/scaffolds_1Mb)
# prune linked sites by scaffolds in parallel
parallel --tmpdir /home/bferrette/sailfish_cl/population.analysis/sailfish/ld_pruning2/tmp \
-j $3 "perl -I /home/bferrette/perl5/lib/perl5 ${ngsld}/scripts/prune_graph.pl --in_file {} --max_kb_dist 50 --min_weight 0.1 --out {}.pruned &> {}.prune_graph.log" ::: $2/snps.ld.*
# concatenate and sort unlinked sites into a single file
cat $(ls $2/snps.ld.*.pruned) | sort -V | sed 's/:/_/' > $2/snps.ld.pruned # input ANGSD
sed -r 's/(.*)_/\1\t/g' snps.ld.pruned > ld_pruned.sites
angsd sites index ld_pruned.sites
# concatenate 'prune_graph.pl' log files
cat $(ls -v $2/snps.ld.*.prune_graph.log) > $2/prune_graph.log
# remove intermediate files
# wc -l $2/snps.ld.scaffold*.pruned
# wc -l $2/snps.ld.pruned
rm $2/snps.ld.HiC_scaffold_*
pigz --best snps.ld
# run stop!!!!

# set output file name
#pruned_gl=$2/snps.ld_pruned.beagle
# extract SNPs matching the '*.ld.pruned' file from the GL file
#gl=$(find $1 -name '*.beagle.gz')
#gl=$(find /home/bferrette/sailfish_denovo/assembly_denovo/snp_calling -name '*.beagle.gz')
#zcat ${gl} | grep -F -f $2/snps.ld.pruned > ${pruned_gl}.tmp
# get linked SNPs inadvertently extracted due to partial string matching
#pigz $2/snps.ld &
#awk '{ print $1 }' ${pruned_gl}.tmp | diff $2/snps.ld.pruned - | grep -Eo 'scaffold.+' > $2/sites2remove
# wc -l $2/sites2remove
# extract the header of the original GL file
#zcat ${gl} | head -1 > ${pruned_gl}
# extract only the unlinked SNPs from the GL (.tmp) file
#grep -v -F -f $2/sites2remove ${pruned_gl}.tmp >> ${pruned_gl}
#awk '{ print $1 }' ${pruned_gl} | diff $2/snps.ld.pruned - | grep -Eo 'scaffold.+'
#grep '' ${pruned_gl}.tmp >> ${pruned_gl}
#sed -i 's/marker/scaffold0_marker/' ${pruned_gl}
#sort -V -k 1,1 ${pruned_gl} | sed 's/scaffold0_marker/marker/' | gzip > ${pruned_gl}.gz
# remove intermediate files
#rm ${pruned_gl} ${pruned_gl}.tmp sites2remove

#for pattern in site1 site2 siteN; do
#  grep "${pattern}\t" ${pruned_gl}.tmp >> ${pruned_gl}
#done

