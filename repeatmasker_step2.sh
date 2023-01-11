# Variables
asm=~/sailfish/repeatmasker_step1/sailfish_assembly_genbank.fasta.masked
dir=~/sailfish/repeatmasker_step2

# RepestMasker
RepeatMasker -s -engine 'ncbi' -pa 32 -lib ../repeatmodeler_results/sailfish_assembly.masked.db-families.fa -dir ${dir} -gff ${asm} > ${dir}/repeatmasker_step2.log

# Concatenate outputs from RepeatMasker
tail -n +4 ./sailfish_assembly_genbank.fasta.masked.out | cat ../repeatmasker_step1/sailfish_assembly_genbank.fasta.out - > repeatmasker.combined.out
cut -f 1,2 ../assembly/sailfish_assembly_genbank.fasta.fai > sailfish_assembly_genbank.tsv
~/anaconda3/envs/RMenv/share/RepeatMasker/util/buildSummary.pl -species 'actinopterygii' -genome sailfish_assembly_genbank.tsv -useAbsoluteGenomeSize repeatmasker.combined.out > repeatmasker.summary.tbl

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
awk '$2>=1000000 { print $1 }' sailfish_assembly_genbank.tsv > scaffolds_1Mb
while read scaff; do
    grep -P "$scaff\t" repeats.bed >> repeats_1Mb.scaffs.bed
done < scaffolds_1Mb
awk '$2>=1000000 { print $0 }' sailfish_assembly_genbank.tsv > scaffolds_sizes

# Make a masked fasta for repeatitive regions
bedtools maskfasta -fi <input FASTA> -bed repeats.bed -fo <output FASTA>

seqtk subseq -l 80 sailfish_cl_assembly_masked.fasta scaffolds_1Mb > sailfish_cl_assembly_masked_scaffs1mbp.fasta

grep -c '>' sailfish_assembly_genbank_masked_scaffs1mbp.fasta

# Get non-repetitive regions
bedtools complement -i repeats_1Mb.scaffs.bed -g scaffolds_sizes > no_repeats.bed
