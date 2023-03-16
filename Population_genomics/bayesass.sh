# Thin VCF file
vcftools --vcf sailfish.pruned.vcf --thin 50000 --remove-indels --recode --recode-INFO-all --out sailfish.pruned.thin.vcf

# Convert VCF in Stacks structure file
~/software/stacks-2.60/populations -V sailfish.pruned.thin.vcf.recode.vcf -M popmap -O . -t 24 --structure

# convert stacks structure into BayesAss input file
~/software/file_converters/stacksStr2immanc.pl -o sailfish.immanc -p WCA,SWA,ECA,SWI,SEI,WCP -s sailfish.pruned.thin.vcf.recode.p.structure

# run the atotune BayesAss to find suitable parameters
python3 /home/bferrette/software/BA3-SNPS-autotune-2.1.2/BA3-SNPS-autotune.py -i sailfish.immanc -l 12034 -b 100000 -g 1000000 -r 10
