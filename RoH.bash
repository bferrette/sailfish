# Build a bcf to run 
gzip -cd ../snp_calling/sailfish.beagle.gz | cut -f1 | tail -n +2 | sed 's/_/\t/2' > sites.txt
wc -l sites.txt
head sites.txt

# Run the script genotype_calling.sh

# Detecting runs of homozygosity (RoH)
https://samtools.github.io/bcftools/bcftools.html#roh
http://samtools.github.io/bcftools/howtos/roh-calling.html
https://samtools.github.io/bcftools/howtos/query.html

# See header
bcftools view -h sailfish.bcf | head
bcftools view -h sailfish.bcf | tail

# View the VCF without AF
bcftools view sailfish.bcf | tail -5

# Now add the AF
bcftools +fill-tags sailfish.bcf  -- -t AF | tail -5

# Grep INFO

bcftools view -h sailfish.bcf | grep "INFO"

# Update the header INFO with AF
bcftools view -h sailfish.bcf > hdr.txt

bcftools reheader -h hdr.txt -o sailfish.reheader.bcf sailfish.bcf --threads 4

# create a allele frequencies files from a tab-delimited file containing the columns: CHROM\tPOS\tREF,ALT\tAF.

bcftools query -f "%CHROM\t%POS\t%REF,%ALT\t%AN\t%AC{0}\n" ../snps.filtered.bcf | awk 'OFS="\t" { print $1,$2,$3,$5/$4 }' | bgzip -c > freqs.tab.gz
tabix -s1 -b2 -e2 freqs.tab.gz

# bcftools RoH all samples
for sample in SFA1-LAT SFA2-LAT SFA3-LAT SFA4-LAT SFA5-LAT SFA6-LAT SFA7-LAT SFA8-LAT SFA9-LAT SFA10-LAT SFA1-LGL SFA2-LGL SFA3-LGL SFA4-LGL SFA5-LGL SFA6-LGL SFA7-LGL SFA8-LGL SFA9-LGL SFA10-LGL; do
	bcftools roh --estimate-AF "PL,-" -s ${sample} snps.filtered.bcf > ${sample}.roh.txt
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# bcftools RoH cluster 1
for sample in SFA1-LAT SFA2-LAT SFA3-LAT SFA4-LAT SFA5-LAT SFA6-LAT SFA7-LAT SFA8-LAT SFA9-LAT SFA10-LAT SFA3-LGL SFA4-LGL SFA5-LGL SFA6-LGL SFA7-LGL; do
	bcftools roh --estimate-AF "PL,cluster1.list" -s ${sample} snps.filtered.bcf > ${sample}.roh.txt
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# bcftools RoH cluster 2
for sample in SFA1-LGL SFA2-LGL; do
	bcftools roh --estimate-AF "PL,cluster2.list" -s ${sample} snps.filtered.bcf > ${sample}.roh.txt
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# bcftools RoH cluster 3
for sample in SFA8-LGL SFA9-LGL SFA10-LGL; do
	bcftools roh --estimate-AF "PL,cluster3.list" -s ${sample} snps.filtered.bcf > ${sample}.roh.txt
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# sailfish populational high
# bcftools RoH cluster atlantic per sample
for sample in SFA01 SFA02 SFA03 SFA04 SFA05 SFA06 SFA07 SFA08 SFA09 SFA10 SFA11 SFA12 SFA13 SFA14 SFA15 SFA16 SFA17 SFA19 SFA20 SFA21 SFA22 SFA23 SFA24 SFA25 SFA26 SFA27 SFA28 SFA29 SFA30 SFA31 SFA32; do
	bcftools roh --estimate-AF "PL,ATL.cluster" -s ${sample} snps.filtered.bcf > ${sample}.roh.txt
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# bcftools RoH cluster indo-pacific per sample
for sample in SFA33 SFA34 SFA35 SFA36 SFA37 SFA38 SFA39 SFA40 SFA41 SFA42 SFA43 SFA44 SFA45 SFA46 SFA47 SFA48 SFA49 SFA50 SFA51 SFA52 SFA53 SFA54 SFA55 SFA56 SFA57 SFA58 SFA59 SFA60 SHS-SA01; do
	bcftools roh --estimate-AF "PL,IDWP.cluster" -s ${sample} snps.filtered.bcf > ${sample}.roh.txt
done
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ${sample}.roh.tsv; done

# bcftools RoH cluster atlantic
bcftools roh --estimate-AF "PL,ATL.cluster" -S ATL.cluster snps.filtered.bcf > ATL.roh.txt
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > ATL.roh.tsv; done

# bcftools RoH cluster indo-pacific
bcftools roh --estimate-AF "PL,IDWP.cluster" -S IDWP.cluster snps.filtered.bcf > IDWP.roh.txt
for f in *.txt; do sample=$(basename $f .roh.txt); grep 'RG' $f | tail -n +2 > IDWP.roh.tsv; done
