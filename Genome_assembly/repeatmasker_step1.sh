#!bin/bash
#
# RepeatMasker v4.1.2
# RepeatMasker uses a sequence search engine to perform it's search for repeats
# https://www.repeatmasker.org/RepeatMasker/
#
# Variables
asm=${HOME}/sailfish_cl/assembly_cl/sailfish_CL_gap_closed_2.fasta
dir=${HOME}/bferrette/sailfish_cl/assembly_cl/repetitive_regions/repeatmasker_step1

# RepestMasker
${HOME}/software/RepeatMasker-4.1.2/RepeatMasker -s -engine 'ncbi' -pa 16 -species 'actinopterygii' -dir ${dir} -gff ${asm} > ${dir}/repeatmasker_step1.log

