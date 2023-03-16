#!/bin/bash
#
# Scaffolding with Hi-C or Omni-C
#
# Requirements:
#  BWA - https://github.com/lh3/bwa
#  Juicer: A One-Click System for Analyzing Loop-Resolution Hi-C Experiments (https://github.com/aidenlab/juicer)
#  3D-DNA: 3D de novo assembly (3D DNA) pipeline (https://github.com/aidenlab/3d-dna)
#  TGS-GapCloser: A gap-closing software tool that uses long reads to enhance genome assembly (https://github.com/BGI-Qingdao/TGS-GapCloser)
#  seqtk: Toolkit for processing sequences in FASTA/Q formats (https://github.com/lh3/seqtk)

# Run the pipeline:
  #Step 1a: copy the assembly you want to scaffold into the designated juicer "reference" dirctory. Please check the juicer manual for help
  #Step 2b: index the assembly with bwa index (see Step 5a)
  #Step 3c: run the generate_site_positions.py script on the assembly in case you use Hi-C data. the restriction enzyme will very likely be MboI; # in case you use Omni-C I would still run it as if you use omni-c to be able to generate the chromosome sizes in the next step
  #Step 4d: after running generate_site_positions, run the follwing command to calculate chromosome sizes:
    awk 'BEGIN{OFS="\t"}{print $1, $NF}' mygenome_myenzyme.txt > mygenome.chrom.sizes
  #Step 5e: run juicer.sh #please hard code the reference and the restriction enzyme into your script. setting parameters does not always work with this script.
    juicer.sh -z <reference> -p <calculated.chrom.sizes> -y <restriction sites, use only for Hi-C> -t <threads> -e    # -e for early exit, we do not need to let it run through completely
  #Step 6f: run 3D-DNA on the merged.nodups.txt file from juicer:
    run-asm-pipeline.sh <reference> <merged_nodups.txt>
  #Step 7g: scaffolding QC:
    #Step 7g1: Quast see Step 7a
    #Step 7g2: check contact-map (.hic) with juicebox --> do you see the correct number of chromosomes?
    #Step 7g3: BUSCO see Step 7b

#Step 8: Gap-closing: to fill gaps (Ns) in your scaffolds with actual sequencing you can use TGS-GapCloser. This can be done multiple times until you see no improvements. The script includes racon and pilon polishing steps to polish the raw reads before gap-closing. long reads must be converted to fasta format!
  #Step8a: convert long-read fastq to fasta:
    seqtk seq -a IN.fastq > OUT.fasta
  #Step8b: run TGS-GapCloser without error correction
    ${HOME}/software/TGS-GapCloser/TGS-GapCloser.sh \
      	--scaff ${HOME}/sailfish/assembly/Istiophprus_platypetrus.fasta \
      	--reads ${HOME}/billfishes/sailfish_cl/data/pacbio/istiophorus_pacbio_clr.fasta \
      	--output sailfish_gapclosed \
      	--ne \
      	--racon ${HOME}/software/racon/build/bin/racon \
      	--r_round 3 \
      	--pilon ${HOME}/software/pilon/pilon-1.23.jar \
      	--p_round 3 \
      	--minmap_arg '-x map-pb -H' \
      	--ngs ${HOME}/billfishes/sailfish_cl/processing/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_1.fq.gz \
      	--ngs ${HOME}/billfishes/sailfish_cl/processing/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_2.fq.gz \
      	--samtools ${HOME}/software/samtools-1.16/samtools \
      	--java ${HOME}/software/jdk-11.0.14/bin/java \
      	--tgstype pb \
      	--thread 16 \
      	--pilon_mem 256G \
      	> TGS-GapCloser.log 2> TGS-GapCloser.err
