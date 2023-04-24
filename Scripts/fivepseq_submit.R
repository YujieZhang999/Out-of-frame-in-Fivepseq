
# This code including all steps to analyse the out-of-frame data in the manuscript

## Step 1: setup the folder
workdir <- "out-of-frame/"
codeDir <- file.path(workdir,"code")
source_code <- file.path(codeDir,"source_code")
fivepseq_count <- file.path(workdir,"fivepseq")
output <- file.path(workdir,"output_data")
rds <- file.path(workdir,"rds")
plotdir <- file.path(workdir,"plot")
whichsample <- "CSM_ctrl"

## Step2: read the count for each gene at each position
### Input data is the output from fivepseq package, output data is saved as FrameCount_List_df for each sample
source(file.path(source_code,"ReadFivepseq_Frame_individual.R"))

## Step3: merge replicates 
### Plot density plot compared with YPD
source(file.path(source_code,"CalculateFSindex.R"))
