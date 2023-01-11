#Merqury
# Activate conda environment
conda activate MERQURYenv

# Best K-mer number
sh ~/anaconda3/envs/MERQURYenv/share/merqury/best_k.sh 618640819
genome: 618640819
tolerable collision rate: 0.001
19.5844

k=19.5844
for sample in SFA7-LAT_1.fastq.gz SFA7-LAT_2.fastq.gz
do
# 1. Build meryl dbs
    meryl k=$k count output $sample.meryl $sample threads=24 memory=48g
done

# 2. Merge
meryl union-sum output sailfish.meryl SFA7-LAT*.meryl

# Run merqury
merqury.sh sailfish.meryl ../sailfish_assembly_genbank.fasta sailfish.merqury.out

# Chromosome-level merqury
# Activate conda environment
conda activate MERQURYenv

# Best K-mer number
sh ~/anaconda3/envs/MERQURYenv/share/merqury/best_k.sh 619036510
genome: 619036510
tolerable collision rate: 0.001
19.5849

k=19.5849
for sample in /home/bferrette/sailfish_cl/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_1.fq.gz /home/bferrette/sailfish_cl/trimmed/SFA29_FDSW202536543-1r_HKFGWDSXY_L4_2.fq.gz
do
# 1. Build meryl dbs
    meryl k=$k count output $sample.meryl $sample threads=24 memory=48g
done

# 2. Merge
meryl union-sum output sailfish.meryl ./SFA29_FDSW202536543-1r_HKFGWDSXY_L4_*.meryl

# Run merqury
merqury.sh ./sailfish.meryl /home/bferrette/sailfish_cl/assembly_cl/sailfish_CL_gap_closed_2.fasta sailfish.merqury.out
