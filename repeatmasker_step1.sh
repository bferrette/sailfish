# Variables
asm=~/sailfish/assembly/sailfish_assembly_genbank.fasta
dir=~/sailfish/repeatmasker_step1

# activate conda environment
conda activate RMenv

# RepestMasker
RepeatMasker -s -engine 'ncbi' -pa 32 -species 'actinopterygii' -dir ${dir} -gff ${asm} > ${dir}/repeatmasker_step1.log
