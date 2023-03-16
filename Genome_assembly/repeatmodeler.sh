#!/bin/bah
#
# RepeaModeler v2.0.2a
# RepeatModeler is a de novo transposable element (TE) family identification and modeling package
# https://www.repeatmasker.org/RepeatModeler/

#Variable
masked_asm=${HOME}/sailfish_cl/assembly_cl/repetitive_regions/repeatmasker_step1/sailfish_CL_gap_closed_2.fasta.masked

#Build database
${HOME}/software/RepeatModeler-2.0.2a/BuildDatabase -engine 'ncbi' -name sailfish.masked.db ${masked_asm}

#RepeatModeler
${HOME}/software/RepeatModeler-2.0.2a/RepeatModeler -engine 'ncbi' -pa 16 -database sailfish.masked.db -LTRStruct > repeatmodeler.log

