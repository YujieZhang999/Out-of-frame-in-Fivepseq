# Out-of-frame analysis in Fivepseq

## Data preprocessing
* using fivepseq package v1.2.0 https://github.com/lilit-nersisyan/fivepseq
* submit ```Preprocess_Fivepseq_yeast.sh``` for fivepseq in yeast for the pipline including cutadapter, extract UMI, align to genome and deduplicate.
* submit ```Preprocess_Human_riboseq.sh```for riboseq in human(Downloaded from GEO: GSE113751, DOI:10.1016/j.molcel.2018.06.041)
* Input data is fastq , output data directly from fivepseq package 
* use the output for Fig1


## 1. Read the count for each gene at each position
```
source(file.path(source_code,"ReadFivepseq_Frame_individual.R"))
```
* Input data is the output from fivepseq package
* output data is saved as FrameCount_List_df for each sample with the frame counts for each gene


## 2. Calculate the gene frameshift index
```
source(file.path(source_code,"CalculateGene_FSindex.R"))
```
* Input data is the FrameCount_List_df 
* use the output for Fig2


## 3. Calculate RNA half-life using SLAM-seq data
* Input data is the T-C conversion counts txt from ```Preprocessing_SLAM-Seq.sh```
* use the ```SLAMSeq_halflife.R``` to calculate the RNA half-life
* use the output for Fig3


## 4. Calculate the codon frameshift index
```
source(file.path(source_code,"CalculateCodon_FSindex.R"))
```
* Input data is codon_pauses.txt (from fivepseq package output) to calculate the frameshift for each codon
* use the output for Fig4 



## 5. Protein abundance change by MS 
* Input data is table for protein quantification statistics
* use ```MS_protein_abundance.R``` to compare the protein change in wt and upf1 mutant
* use the output for Fig5



## 6. Count the out-of-frame genes in bacteria and human
* Input data is the frameshift index calculated from ```CalculateGene_FSindex.R```
* use ```Out_of_frame_bacteria_human.R``` to correlate the frameshift index with codon composition
* use the output for Fig6

