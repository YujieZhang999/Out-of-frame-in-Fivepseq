#!/bin/sh
#SBATCH -A 
#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 2-50:00:00
#SBATCH -J riboseq_preprocess_and_fivepseq
#SBATCH -o log/riboseq_preprocess_and_fivepseq.log


genomedir="/crex/proj/snic2019-30-56/lilit/genome/human/ensembl"

/crex/proj/snic2019-30-56/nobackup/projects/Yujie_HT5Pseq/Ribosome_profiling/2018_AA_human/src/prefivepseq_small-RNA_modify.sh \
    -f fastq \
    -o preprocess_riboseq \
    -a $genomedir/gff/Homo_sapiens.GRCh38.93.gff3 \
    -s du \
    -g $genomedir/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa \


echo "########################################################"
echo "####   fivepseq count_and_plot on yeast "
echo "########################################################"

fivepseq \
    -g $genomedir/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -a $genomedir/gff/Homo_sapiens.GRCh38.93.gff3 \
    --conflicts add \
    -b "preprocess_riboseq/star_align_subset/*.bam" \
    -o fivepseq_riboseq_subset \
    -gf "human_high_fs_genes.txt"
    -t riboseq \
    -queue 30 \
    --oof
