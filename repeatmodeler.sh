#Variable
masked_asm=../repeatmasker_step1/sailfish_assembly_genbank.fasta.masked
#Build database
BuildDatabase -engine 'ncbi' -name sailfish_assembly.masked.db ${masked_asm}
#RepeatModeler
~/software/RepeatModeler-2.0.1/RepeatModeler -engine 'ncbi' -pa 32 -database sailfish_assembly.masked.db -LTRStruct > repeatmodeler.log
