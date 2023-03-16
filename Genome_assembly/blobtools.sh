#!/bin/bash
#
# Blobtools: Modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets
# https://github.com/DRL/blobtools
#
# Requirements:
#  Blobtools v1.1
#  BWA - https://github.com/lh3/bwa
#  Minimap2 - https://github.com/lh3/minimap2
#  NCBI-BLAST - https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
#
# map reads to your reference (bwa for short-reads ; minimap2 for long reads (use -a option to generate SAM instead of paf)):
    minimap2 -ax map-pb -H ../Istiophorus_platypterus.fasta ../sailfish_pacbio.fasta -t 12 | samtools sort -l 9 -@ 16 -o sailfish.pacbio.bam -O BAM && samtools index sailfish.pacbio.bam
  #Step 1d: Blobtools to check for contigs with contamination
    #Step 1d1: we will need BAM-files for each read type --> you can use the BAM-files from Step 7c
    #Step 1d2: in addition we need to blast our assembly against the nucleotide database from ncbi (update the nt database first if necessary. This can take a long time):
      ${HOME}/software/ncbi-blast-2.13.0+/bin/blastn -task megablast -query ../sailfish_CL_gap_closed_2.fasta -db /opt/software/blastdb/nt/nt -outfmt '6 qseqid staxids bitscore std' -num_threads 16 -evalue 1e-25 -out sailfish_cl.blast.hits
    #Step 1d3: run blobtools, I used blobtools v1.1 so the command can vary with a newer version:
      #Step 1d3.1: create BlobDB:
      # wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
      # mkdir data && mv taxdump.tar.gz data/
      # tar -zxf taxdump.tar.gz -C ./ nodes.dmp names.dmp
      blobtools create \
      -i ../Istiophorus_platypterus.fasta \
      -b ../minimap2/sailfish.pacbio.bam \
      -b ${HOME}/sailfish/processing/mapping.shortreads/SFA29_FDSW202536543-1r_HKFGWDSXY_L4.sorted.bam \
      -t ../blastn/sailfish_cl.blast.hits \
      -o blobtools \
      --nodes data/nodes.dmp \
      --names data/names.dmp \
      &> blobtools.log
      #Step 1d3.2: create table from BlobDB to find contaminated contigs:
        blobtools view -i sailfish.blobtools.blobDB.json -o sailfish_cl.blobtools
      #Step 1d3.3: generate plots:
        blobtools plot -i sailfish.blobtools.blobDB.json -o sailfish_cl.blobtools.plot --rank order --format svg

