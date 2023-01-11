#!/bin/bash
# Description: Run the ADMIXTURE software http://dalexander.github.io/admixture/
# create input file https://speciationgenomics.github.io/ADMIXTURE/
# count the number contigs (optional) - awk '{ print $1 }' ~/sailfish/sailfish_assembly_racon_3_pilon_3.fasta.fai | wc -l
plink1.9 --bcf sailfish.filtered.bcf --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --maf 0.025 --hwe 1e-6 midp --indep-pairwise 50 10 0.1 --make-bed --out plink.filtered
plink1.9 --bfile plink.filtered --extract plink.filtered.prune.in --keep-allele-order --allow-extra-chr --make-bed --out sailfish.LDprune
plink1.9 --bfile sailfish.LDprune --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --out sailfish.pruned

# ADMIXTURE does not accept chromosome names that are not human chromosomes.
# We will thus just exchange the first column by 0
awk '{$1=0;print $0}' sailfish.LDprune.bim > sailfish.LDprune.bim.tmp
mv sailfish.LDprune.bim.tmp sailfish.LDprune.bim

# Run ADMIXTURE for different values of K and compile results on cross-validation error across values
prefix=sailfish.LDprune
Klow=1
Khigh=10
for ((K=$Klow;K<=$Khigh;K++)); do admixture --cv=1000 $prefix.bed $K -j48| tee log.$prefix.$K.out; done
echo '# CV results' > $prefix.CV.txt
for ((K=$Klow;K<=$Khigh;K++)); do awk -v K=$K '$1=="CV"{print K,$4}' log.$prefix.$K.out >> $prefix.CV.txt; done

# Make a plot of the cross-validation error as a function of K in R (needs improvement of graph)
# http://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization
tbl <- read.table('sailfish.LDprune.CV.txt')
par(mar = c(4,4,2,2), cex=0.7)
plot(tbl$V1, tbl$V2, xlab='K', ylab='Cross-validation Error', pch=16, type='l')

# Run ADMIXTURE for different values of K, each with 10 replicates
prefix=sailfish.LDprune
for r in {1..10}; do for K in {1..6};
do admixture --seed=RANDOM ${prefix}.bed $K -j24 -B2000
mv ${prefix}.${K}.Q ${prefix}.K${K}r${r}.Q
done; done
