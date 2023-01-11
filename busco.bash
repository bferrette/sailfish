# BUSCO sailfish

# https://gitlab.com/ezlab/busco/-/tree/master

# Run BUSCO
busco -m genome -i ~/sailfish/sailfish_assembly_racon_3_pilon_3.fasta -o busco -l actinopterygii_odb10 --cpu 40

# BUSCO plot generation tool.
# Place all BUSCO short summary files (short_summary.[generic|specific].dataset.label.txt) in a single folder.
# It will be your working directory, in which the generated plot files will be written.

# usage: python3 generate_plot.py -wd [WORKING_DIRECTORY] [OTHER OPTIONS]

python3 ./scripts/generate_plot.py -wd ./ -rt specific
