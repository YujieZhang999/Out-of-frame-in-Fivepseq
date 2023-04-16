

#!/bin/bash -l

f_dir='fastq'
fname=${f_dir##*/}

echo ""
echo "$f_dir"
echo "$fname"

export projectdir=/project/Fivepseq
export refdir=$projectdir/genome/SC


/proj/snic2019-30-56/lilit/prefivepseq/prefivepseq_yeast.sh \
    -f $f_dir \
    -i $refdir/star_index \
    -o preprocess/$fname \
    -a $refdir/gff/* \
    -g $refdir/fa/* \
    > log/$fname.log \



fivepseq \
    -g $refdir/fa/* \
    -a refdir/gff/* \
    --conflicts add \
    -b ../preprocess/fastq/align_dedup/'*R1_001.fastqAligned.primary.out_dedup.bam' \
    -o ../fivepseq/ \
    -t FS \
    -queue 30 \
    --oof
