#Quast
conda activate QUASTenv
# Run
quast.py \
        -o ./quast_results \
        --min-contig 1 \
        --thread 12 \
        --split-scaffolds \
        --large \
        --contig-thresholds 1000,10000,100000,1000000,10000000 \
        Istiophorus_platypterus_ASM1685934v1.fasta \
        sailfish_CL_gap_closed_2.fasta
