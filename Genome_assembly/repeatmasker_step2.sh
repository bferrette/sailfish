#!bin/bash
#
# RepeatMasker v4.1.2
# RepeatMasker uses a sequence search engine to perform it's search for repeats
# https://www.repeatmasker.org/RepeatMasker/
#
# Variables
asm=${HOME}/sailfish_cl/assembly_cl/repetitive_regions/repeatmasker_step1/sailfish_CL_gap_closed_2.fasta.masked
dir=${HOME}/sailfish_cl/assembly_cl/repetitive_regions/repeatmasker_step2

# RepestMasker
${HOME}/software/RepeatMasker-4.1.2/RepeatMasker -s -engine 'ncbi' -pa 16 -lib ../repeatmodeler/sailfish.masked.db-families.fa -dir ${dir} -gff ${asm} > ${dir}/repeatmasker_step2.log

# Concatenate outputs from RepeatMasker
tail -n +4 ./sailfish_CL_gap_closed_2.fasta.masked.out | cat ../repeatmasker_step1/sailfish_CL_gap_closed_2.fasta.out - > repeatmasker.combined.out
cut -f 1,2 ../../../assembly_cl/sailfish_CL_gap_closed_2.fasta.fai > sailfish_assembly_cl.tsv
${HOME}/software/RepeatMasker-4.1.2/util/buildSummary.pl -species 'actinopterygii' -genome sailfish_assembly_cl.tsv -useAbsoluteGenomeSize repeatmasker.combined.out > repeatmasker.summary.tbl

# Create a BED file of identified repeated regions while merging overlapping and adjacent repepats
tail -n +4 repeatmasker.combined.out \
| sed 's/^\s*//' \
| sed -E 's/\s+/\t/g' \
| cut -f 5-7 \
| awk 'OFS="\t" { print $1, $2-1, $3 }' \
| sort -V \
| bedtools merge \
> ./repeats.bed

# Remove scaffolds smaller than 1Mb
awk '$2>=1000000 { print $1 }' sailfish_assembly_cl.tsv > scaffolds_1Mb
while read scaff; do
    grep -P "$scaff\t" repeats.bed >> repeats_1Mb.scaffs.bed
done < scaffolds_1Mb
awk '$2>=1000000 { print $0 }' sailfish_assembly_cl.tsv > scaffolds_sizes

# Make a masked fasta for repeatitive regions
bedtools maskfasta -fi ${HOME}/sailfish_cl/assembly_cl/sailfish_CL_gap_closed_2.fasta -bed repeats.bed -fo sailfish_cl_masked.fasta

seqtk subseq -l 80 sailfish_cl_masked.fasta scaffolds_1Mb > sailfish_cl_masked_scaffs1mbp.fasta

grep -c '>' sailfish_cl_masked_scaffs1mbp.fasta

# Get non-repetitive regions
bedtools complement -i repeats_1Mb.scaffs.bed -g scaffolds_sizes > no_repeats.bed

