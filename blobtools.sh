# Blobtools - 
# map reads to your reference (bwa for short-reads -see Step 5a-c-; minimap2 for long reads (use -a option to generate SAM instead of paf)):
    minimap2 -ax map-pb ../sailfish_CL_gap_closed_2.fasta ../sailfish_pacbio.fasta -t 12 | samtools sort -@ 12 -o sailfish.pacbio.bam -O BAM && samtools index sailfish.pacbio.bam
  #Step 7d: Blobtools to check for contigs with contamination
    #Step 7d1: we will need BAM-files for each read type --> you can use the BAM-files from Step 7c
    #Step 7d2: in addition we need to blast our assembly against the nucleotide database from ncbi (update the nt database first if necessary. This can take a long time):
      /home/bferrette/software/ncbi-blast-2.13.0+/bin/blastn -task megablast -query ../sailfish_CL_gap_closed_2.fasta -db /opt/software/blastdb/nt/nt -outfmt '6 qseqid staxids bitscore std' -num_threads 12 -evalue 1e-25 -out sailfish_cl.blast.hits
    #Step 7d3: run blobtools, I used blobtools v1.1 so the command can vary with a newer version:
      #Step 7d3.1: create BlobDB:
      # wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
      # mkdir data && mv taxdump.tar.gz data/
      # tar -zxf taxdump.tar.gz -C ./ nodes.dmp names.dmp
      blobtools create \
      -i ../sailfish_CL_gap_closed_2.fasta \
      -b ../minimap2/sailfish.pacbio.bam \
      -b /home/bferrette/sailfish_cl/processing/mapping.shortreads/SFA29_FDSW202536543-1r_HKFGWDSXY_L4.sorted.bam \
      -t ../blastn/sailfish_cl.blast.hits \
      -o sailfish_cl.blobtools \
      --nodes data/nodes.dmp \
      --names data/names.dmp \
      &> blobtools.log
      #Step 7d3.2: create table from BlobDB to find contaminated contigs:
        blobtools view -i sailfish_cl.blobtools.blobDB.json -o sailfish_cl.blobtools
      #Step 7d3.3: generate plots:
        blobtools plot -i sailfish_cl.blobtools.blobDB.json -o sailfish_cl.blobtools.plot --rank order --format svg
