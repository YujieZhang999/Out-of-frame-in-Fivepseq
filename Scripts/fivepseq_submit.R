
# This code includes the analysis of out-of-frame data in the manuscript

## Setup the folder
workdir <- "out-of-frame/fivepseq"
codeDir <- file.path(workdir,"code")
source_code <- file.path(codeDir,"source_code")
fivepseq_count <- file.path(workdir,"fivepseq_count")
output <- file.path(workdir,"output_data")
rds <- file.path(workdir,"rds")
plotdir <- file.path(workdir,"plot")
input <- file.path(workdir,"input_data")
whichsample <- "CSM_ctrl"

## Read the count for each gene at each position
### Input data is the output from fivepseq package
### output data is saved as FrameCount_List_df for each sample
source(file.path(source_code,"ReadFivepseq_Frame_individual.R"))

## Calculate the gene frameshift index
### Input data is the FrameCount_List_df
### use the output for Fig2
source(file.path(source_code,"CalculateGene_FSindex.R"))

## Calculate RNA half-life using SLAM-seq data
### Input data is the T-C conversion counts txt from Preprocessing_SLAM-Seq.sh
### use the SLAMSeq_halflife.R to calculate the half-life
### use the output for Fig3


## Calculate the codon frameshift index
### Input data is codon_pauses.txt from the output from fivepseq package
### use the output for Fig4
source(file.path(source_code,"CalculateCodon_FSindex.R"))


## Count the out-of-frame genes in bacteria and human
### Using Out_of_frame_bacteria_human.R to 