#!/bin/bah
#
# Repeat Landscapes - https://littlebioinformatician.wordpress.com/scripts/repeat-landscape/
perl ${HOME}/software/RepeatMasker-4.1.2/util/calcDivergenceFromAlign.pl -s sailfish.divsum ${HOME}/sailfish_cl/assembly_cl/repetitive_regions/repeatmasker_step2/sailfish_CL_gap_closed_2.fasta.masked.cat.gz
perl ${HOME}/software/RepeatMasker-4.1.2/util/createRepeatLandscape.pl -g 612904761 -div sailfish.divsum > sailfish.landscape.html

