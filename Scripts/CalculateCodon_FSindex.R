#Setup
library("tidyverse")
library(dplyr)
library("zoo") 
library("ggpubr")
opts <- options(stringsAsFactors = F)

# # -----------------------------------------------------------------------
# Input data is codon_pauses.txt from the output from fivepseq package
# # -----------------------------------------------------------------------
codon_file <- function(whichframe="codon_pauses_f0.txt", whichcodon="f0_codon"){
  codon.files <- list.files(fivepseq_count, pattern = "^codon_pauses*",full.names = T, recursive = T)
  codon.files_list <- codon.files[which(basename(codon.files)==whichframe)]
  for (codon.f in codon.files_list) {
    codon <- read.table(codon.f, sep = "\t", header = T,  row.names=1)
    codon <- codon[,c(1:33)]
    colnames(codon) <- seq(-30,2)
    name <- basename(dirname(codon.f))
    mycodon[[name]][[whichcodon]] <- codon
  }
  return(mycodon)
}

codon_file_f0 <- codon_file(whichframe="codon_pauses_f0.txt",whichcodon="f0_codon")
codon_file_f1 <- codon_file(whichframe="codon_pauses.txt",whichcodon="f1_codon")
codon_file_f2 <- codon_file(whichframe="codon_pauses_f2.txt",whichcodon="f2_codon")
codon_list <- list("codon_file_f0"=codon_file_f0,"codon_file_f1"=codon_file_f1,"codon_file_f2"=codon_file_f2)
#saveRDS(codon_list,file.path(rds,paste(whichsample,"codon_frame.list.rds",sep = "_")))

# # -----------------------------------------------------------------------
# Extract codon peaks for each samples and save in the dataframe from list
# # -----------------------------------------------------------------------
GetFrameCount <- function(sample, whichcodon, dataset = "fivepseq_yeast") {
  mysum <- sample[[whichcodon]] %>% as.matrix() %>% sum()
  b <- sample[[whichcodon]] %>% 
    as_tibble() %>%
    dplyr::mutate_all(~ (. / mysum) * 1000000) %>% 
    if (dataset == "fivepseq_yeast") {
      dplyr::select(`-18`, `-17`)  
    } else if (dataset == "fivepseq_bacteria") {
      dplyr::select(`-15`, `-14`)  
    } else if (dataset == "riboseq_yeast") {
      dplyr::select(`-16`, `-15`) 
    } 
  b$codons <- unlist(sapply(rownames(sample[[whichcodon]]), function(x) strsplit(as.character(x), split = "_")[[1]][2]))
  b$aa <- rownames(sample[[whichcodon]])
  # c <- b %>% left_join(tAI, by = "codons") 
  return(b)
}
codon_aver_f1 <- map2(codon_file_f1, "f1_codon", GetFrameCount)

# # -----------------------------------------------------------------------
# Calculate codon frameshift index
# # -----------------------------------------------------------------------
cal_fs <- function(df) {
  val <- log2(df$`-17`/df$`-18`) # for bacteria_fivepseq, use -14( F1) and -15 (F0); 
                                 # for yeast_riboseq , use -15 (F1) and -16 (F0); 
                                 # for yeast_fivepseq, use -17 (F1) and -18 (F0);
  return(val)
  }   

codonFS_df <- map_dfc( .x = codon_aver_f1, .f = ~ cal_fs(.x)) 
codonFS_df <- as.data.frame(codonFS_df) 
rownames(codonFS_df) <- codon_aver_f1[[1]]$aa

samples <- colnames(codonFS_df)
group_condition <- lapply(strsplit(samples, "_"), function(x) (x[1:2])) # it depends on the sample name
group_condition <- sapply(group_condition, paste, collapse = "_") 

# function for sample comparsion
diff_frame <- function(control_ind, sample_ind, FS.df){
  c.df <- na.omit(as.data.frame(FS.df[,control_ind]))
  a.df <- na.omit(as.data.frame(FS.df[,sample_ind]))
  # check if there are replicates or not
  if (length(control_ind) == 1 & length(sample_ind) == 1) {
    diff <- merge(as.data.frame(c.df), as.data.frame(a.df), by = 0)
  } else {
    diff <- merge(as.data.frame(rowMeans(c.df)), as.data.frame(rowMeans(a.df)), by = 0)
  }
  diff <- column_to_rownames(diff, var = "Row.names")
  colnames(diff) <- c("fs_rich","fs_poor")
  diff <- diff %>% mutate("diff" = fs_rich - fs_poor) %>% na.omit()
  return(diff)
}

TwoGroup_Compare <- function(control,sample,FS.df){
  c.ind <- which(grepl(control, samples))
  a.ind <- which(grepl(sample, samples)) 
  diff<- diff_frame(c.ind,a.ind,FS.df)
  return(diff)
}

# Codon frameshift compare and merge the replicates
frameCodon.mat <- matrix(data = 0, nrow = 64, ncol = length(unique(group_condition)))
rownames(frameCodon.mat) <- rownames(codonFS_df)
colnames(frameCodon.mat) <- unique(group_condition)

for (i in 1:length(unique(group_condition))){
  df_Comp <- TwoGroup_Compare("YPD_ctrl",unique(group_condition)[i],codonFS_df) # Choose the control sample
  name <- rownames(df_Comp)
  frameCodon.mat[name,i] <- df_Comp[,"diff"] # Extracted the useful column
}
# save df_Comp and used it for the Fig5.R plot














