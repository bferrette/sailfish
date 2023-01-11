# renomear nomes dos arquivos
for f in *.bam; do echo "$f ${f%_*_*_*.clean.bam}.clean.bam"; done

for f in *.bam; do
  mv $f ${f%_*_*_*.clean.bam}.clean.bam; done

for f in *.bam; do
  echo "$f ${f%.clean.bam}"; done

for f in *.bai; do
  mv $f ${f%_*_*_*.clean.bam.bai}.clean.bam.bai; done

parallel -j 12 "grep -P {}: snps.ld > snps.ld.{}" ::: $(cat /home/bferrette/sailfish/repeatmasker_step2/scaffolds_1Mb)

wget -O- http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib
eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`
echo 'eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`' >> ~/.bash_profile
echo 'export MANPATH=$HOME/perl5/man:$MANPATH' >> ~/.bash_profile

parallel -j 20 "perl -I /home/bferrette/perl5/lib/perl5 ${ngsld}/scripts/prune_graph.pl --in_file {} --max_kb_dist 20 --min_weight 0.1 --out {}.pruned &> {}.prune_graph.log" ::: snps.ld.*

cat $(ls snps.ld.*.pruned) | sort -V | sed 's/:/_/' > snps.ld.pruned
cat $(ls -v snps.ld.*.prune_graph.log) > prune_graph.log

wc -l snps.ld.ctg*.pruned
wc -l snps.ld.pruned
rm snps.ld.scaffold*.pruned
rm snps.ld.scaffold*.log

pruned_gl=snps.ld_pruned.beagle

gl=$(find $1 -name '*.beagle.gz')
gl=$(find /home/bferrette/sailfish_denovo/assembly_denovo/snp_calling -name '*.beagle.gz')
zcat ${gl} | grep -F -f snps.ld.pruned > ${pruned_gl}.tmp
awk '{ print $1 }' ${pruned_gl}.tmp | diff snps.ld.pruned - | grep -Eo 'chr.+' > sites2remove
zcat ${gl} | head -1 > ${pruned_gl}

pigz snps.ld &
awk '{ print $1 }' snps.ld_pruned.beagle.tmp | diff snps.ld.pruned - | grep -Eo 'scaffold.+' > sites2remove
wc -l sites2remove
zcat ${gl} | head -1 > ${pruned_gl}
grep -v -F -f sites2remove ${pruned_gl}.tmp >> ${pruned_gl}
wc -l snps.ld_pruned.beagle

zgrep 'scaffold28_size8315501_1436139' ${gl} >> snps.ld_pruned.beagle
sed -i 's/marker/scaffold0_marker/' snps.ld_pruned.beagle
sort -V -k 1,1 snps.ld_pruned.beagle > snps.ld_pruned.beagle.sorted
sed -i 's/scaffold0_marker/marker/' snps.ld_pruned.beagle
sed -i 's/scaffold0_marker/marker/' snps.ld_pruned.beagle.sorted

rm snps.ld_pruned.beagle
mv snps.ld_pruned.beagle.sorted snps.ld_pruned.beagle
gzip snps.ld_pruned.beagle
rm snps.ld_pruned.beagle.tmp
