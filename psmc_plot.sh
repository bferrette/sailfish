# load file
psmc_istiophorus=$(ls -f {SFA15,SFA60}.psmc | tr '\n' ' ')
psmc_kajikia=$(ls -f WHM01.psmc | tr '\n' ' ')
psmc_makaira=$(ls -f BUM01.merged.psmc | tr '\n' ' ')
psmc_tetrapturus=$(ls -f SPF01.merged.psmc | tr '\n' ' ')
psmc_xiphias=$(ls -f Xiphias_gladius_SAMN16522248_BGISeq.psmc | tr '\n' ' ')

# No bootstrap
psmc_plot.pl -u 2.2e-09 -g 4.3 -G -R -M 'SFA15','SFA60' -P 'right top' -p sailfish.psmc.plot ${psmc_istiophorus}
psmc_plot.pl -u 2.2e-09 -g 5.5 -G -R -M 'WHM01' -P 'right top' -p kajikia.psmc.plot ${psmc_kajikia}
psmc_plot.pl -u 2.2e-09 -g 4.5 -G -R -M 'BUM01' -P 'right top' -p makaira.psmc.plot ${psmc_makaira}
psmc_plot.pl -u 2.2e-09 -g 4.0 -G -R -M 'SPF01.merged' -P 'right top' -p tetrapturus.psmc.plot ${psmc_tetrapturus}
psmc_plot.pl -u 5.0e-08 -g 6.5 -G -R -M 'Xiphias_gladius_SAMN16522248_BGISeq' -P 'right top' -p xiphias.psmc.plot ${psmc_xiphias}

#Booststrap
#GENERATION LENGTH (YEARS): 4.3 years
samples=(SFA01 SFA02 SFA03 SFA04 SFA05 SFA06 SFA07 SFA08 SFA09 SFA10 SFA11 SFA12 SFA13 SFA14 SFA15 SFA16 SFA17 SFA18 SFA19 SFA20 SFA21 SFA22 SFA23 SFA24 SFA25 SFA26 SFA27 SFA28 SFA29 SFA30 SFA31 SFA32 SFA33 SFA34 SFA35 SFA36 SFA37 SFA38 SFA39 SFA40 SFA41 SFA42 SFA43 SFA44 SFA45 SFA46 SFA47 SFA48 SFA49 SFA50 SFA51 SFA52 SFA53 SFA54 SFA55 SFA56 SFA57 SFA58 SFA59 SFA60 SFA61.20x)
for sample in ${samples[@]}; do
  cat ${sample}.psmc ${sample}.round-*.psmc > ${sample}.combined.psmc
  psmc_plot.pl -u 2.2e-09 -g 4.3 -X 0 -R -Y 300 -T "${sample}" -p ${sample}.combined ${sample}.combined.psmc
done

# GENERATION LENGTH (YEARS): 4.5-6.5 years
samples=(WHM01)
for sample in ${samples[@]}; do
  cat ${sample}.psmc ${sample}.round-*.psmc > ${sample}.combined.psmc
  psmc_plot.pl -u 2.2e-09 -g 5.5 -X 0 -R -Y 0 -T "${sample}" -p ${sample}.combined ${sample}.combined.psmc
done

# GENERATION LENGTH (YEARS): 4.5-6 years
samples=(BUM01.merged)
for sample in ${samples[@]}; do
  cat ${sample}.psmc ${sample}.round-*.psmc > ${sample}.combined.psmc
  psmc_plot.pl -u 2.2e-09 -g 4.5 -X 0 -R -Y 0 -T "${sample}" -p ${sample}.combined ${sample}.combined.psmc
done

# GENERATION LENGTH (YEARS): unknown
samples=(SPF01.merged)
for sample in ${samples[@]}; do
  cat ${sample}.psmc ${sample}.round-*.psmc > ${sample}.combined.psmc
  psmc_plot.pl -u 2.2e-09 -g 4 -X 0 -R -Y 0 -T "${sample}" -p ${sample}.combined ${sample}.combined.psmc
done

# GENERATION LENGTH (YEARS): 6.5 years
samples=(Xiphias_gladius_SAMN16522248_BGISeq)
for sample in ${samples[@]}; do
  cat ${sample}.psmc ${sample}.round-*.psmc > ${sample}.combined.psmc
  psmc_plot.pl -u 5.0e-08 -g 6.5 -X 0 -Y 0 -R -G -T "${sample}" -p ${sample}.combined ${sample}.combined.psmc
done
