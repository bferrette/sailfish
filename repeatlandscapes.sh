# repeatlandscapes - https://littlebioinformatician.wordpress.com/scripts/repeat-landscape/

perl /home/bferrette/anaconda3/envs/RMenv/share/RepeatMasker/util/calcDivergenceFromAlign.pl -s sailfish.divsum /home/bferrette/sailfish_cl/assembly_cl/repeatmasker_step2/sailfish_CL_gap_closed_2.fasta.masked.cat.gz

perl /home/bferrette/anaconda3/envs/RMenv/share/RepeatMasker/util/createRepeatLandscape.pl -g 650000000 -div sailfish.divsum > sailfish.landscape.html
