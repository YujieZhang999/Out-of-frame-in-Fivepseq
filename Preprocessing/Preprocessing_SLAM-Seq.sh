# Before start:
# -- Download slamseq pipline from https://nf-co.re/slamseq
# -- prepare slamsample_meta.txt
# -- Reverse Complement reads1 if you are using strand sepecific library preparation

module load bioinfo-tools
module load gnuparallel/20180822
module load SeqKit/0.15.0
module load bioinfo-tools Nextflow/20.10.0

export projectdir=/project/SLAM-seq
export fastqdir=$projectdir/fastq
export workdir=$projectdir/slamdunk
export rcfastqdir=$projectdir/RCread1
export metadir=$projectdir/meta
export utrbed=$projectdir/gtx_update2019_roman.bed
export NXF_HOME=$projectdir/nextflow
export NXF_SINGULARITY_CACHEDIR=$NXF_HOME/singularity



# Step1: Reverse Complement fastq and save to a new folder
cd $fastqdir
find *.fastq.gz | parallel -a - -j 16 seqkit seq --seq-type DNA -r -p {} '|' gzip -c  '>' ../RCread1/{} 

# Step 2: Using nextflow slamdunk perform T to C conversions calling
# for default setting conversion >= 1 is considered as labeled
# output "_filtered_tcount_collapsed.csv" later used for RNA halflife calculation
cd $workdir

nextflow run $NXF_HOME/assets/nf-core/slamseq/main.nf \
-profile uppmax \
--trim5 2 \
--conversions 1 \
--project 'snic2018-8-292' \
--fasta /crex/data/igenomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa \
--bed $utrbed \
--input $metadir/slamsample_meta.txt \
--read_length 150 \
--max_memory '300.GB' \
--max_time '120.h' \
--max_cpus 20 \
--skip_deseq2 \
-bg

