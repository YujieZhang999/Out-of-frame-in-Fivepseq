# Out-of-frame-in-Fivepseq

## Data preprocessing
##### using fivepseq package https://github.com/lilit-nersisyan/fivepseq 
##### submit Preprocess_Fivepseq_yeast.sh for fivepseq in yeast
##### submit Preprocess_Human_riboseq.sh for riboseq in human(GEO: GSE113751)
##### Input data is fastq , output data directly from fivepseq package
##### use the output for Fig1


## Read the count for each gene at each position
##### Input data is the output from fivepseq package
##### output data is saved as FrameCount_List_df for each sample
source(file.path(source_code,"ReadFivepseq_Frame_individual.R"))

## Calculate the gene frameshift index
##### Input data is the FrameCount_List_df
##### use the output for Fig2
source(file.path(source_code,"CalculateGene_FSindex.R"))

## Calculate RNA half-life using SLAM-seq data
##### Input data is the T-C conversion counts txt from Preprocessing_SLAM-Seq.sh
##### use the SLAMSeq_halflife.R to calculate the RNA half-life
##### use the output for Fig3


## Calculate the codon frameshift index
##### Input data is codon_pauses.txt from the output from fivepseq package
##### use the output for Fig4
source(file.path(source_code,"CalculateCodon_FSindex.R"))


## Protein abundance change by MS 
##### Input data is table for protein quantification statistics
##### use MS_protein_abundance.R to compare the protein change in wt and upf1 mutant
##### use the output for Fig5



## Count the out-of-frame genes in bacteria and human
##### Input data is the frameshift index calculated from CalculateGene_FSindex.R
##### use Out_of_frame_bacteria_human.R to correlated the frameshift index with codon composition
##### use the output for Fig6

