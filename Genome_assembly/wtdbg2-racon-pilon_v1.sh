#!/bin/bash
#
# wtdbg2-racon-pilon.pl: Wrapper for overlap assembly of long- and short-reads and polishing
# https://github.com/schellt/wtdbg2-racon-pilon
#
# Run the pipeline
wtdbg2-racon-pilon.pl \
    -l ${HOME}/billfishes/sailfish_cl/data/pacbio/istiophorus_pacbio_clr.fa.gz \
    -x rsII \
    -racon-rounds 3 \
    -p ${HOME}/billfishes/sailfish_cl/processing/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_1.fq.gz,${HOME}/billfishes/sailfish_cl/processing/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_2.fq.gz \
    -pilon-path ${HOME}/software/pilon/pilon-1.23.jar \
    -pilon-rounds 3 \
    -xmx 256g \
    -t 32 \
    -o . \
    -racon-opts '-u' \
    -bwa-opts '-a -c 10000' \
    -wtdbg-opts '-g 700m' \
    -pilon-opts '--diploid' \
    -v
