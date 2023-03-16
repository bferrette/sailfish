#!/bin/bash
#
# create a submission file for each annotated mitogenome
#
# create tbl files for each annoatation file from MitoAnnotator (http://mitofish.aori.u-tokyo.ac.jp/annotation/input.html)
# adapted version of the python script to Python 3 from https://github.com/chrishah/MitoFish2tbl
#
# create a list with the names of each annotated mitogenome
ls -1 *.fa | grep -v 'genes' | cut -f1 -d '.' > mitogenomeslist.txt

# create a 5-column feature table by converting the MitoAnnotator annotation file to NCBI feature table format (*.tbl)
for mitogenome in $(cat mitogenomeslist.txt)
do
python3 ./mitofish2tbl.py ${mitogenome}.txt
done

# create the genbank submission files:
# create a GenBank Submission Template https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
# table2asn is a command-line program that creates sequence records for submission to GenBank. https://www.ncbi.nlm.nih.gov/genbank/table2asn/
# run the command for every annotated mitogenome in the directory
for mitogenome in $(cat mitogenomeslist.txt)
do
/opt/table2asn -i ${mitogenome}.fa -f ${mitogenome}.tbl -t template.sbt -euk -a s -V bv -T -Z -j "[mgcode=2] technique [tech=wgs] [topology=circular] [location=mitochondrion] [organism=Istiophorus platypterus]"
done

# after to check the BioProject and BioSample of each annotated mitogenome,
# directly submittedthe files *.sqn or *.gbf to GenBank by electronic mail at: gb-sub@ncbi.nlm.nih.gov

