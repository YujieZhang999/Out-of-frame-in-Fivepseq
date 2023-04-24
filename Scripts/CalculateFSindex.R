# Input data is FrameCount_List_df from ReadFivepseq_Frame_individual.R script
# Extract f.change from list
# for each sample, save the f.change in the matrix list

FrameCount_List_df <- readRDS(file.path(rds,paste(whichsample,"IndiV_FrameCount_list.rds",sep = "_")))
samples = names(FrameCount_List_df)

# Read transcript_descriptors.txt and check the dimentionality
t.des.name <- "transcript_descriptors.txt"
t.des.files <- list.files(fivepseq_count, t.des.name, full.names = T, recursive = T)
t.des <- read.table(t.des.files[[1]], sep = "\t", header = T)
Nrow <- dim(t.des)[1]

# Define expected values from transcript_descriptor file
Exp_yeast_Nrow <- 6600
Exp_human_Nrow <- 82724
Exp_bsubs_Nrow <- 4325
Exp_lpla_Nrow <- 3063

# Checkpoint test
if (Nrow == Exp_yeast_Nrow) {
  cat("Passed yeast checkpoint test!\n")
} else if (Nrow == Exp_human_Nrow) {
  cat("Passed human checkpoint test!\n")
} else if (Nrow == Exp_bsubs_Nrow) {
  cat("Passed b.subs checkpoint test!\n")
} else if (Nrow == Exp_lpla_Nrow) {
  cat("Passed l.pla checkpoint test!\n")
} else {
  cat("Failed checkpoint test. Nrow value does not match expected values.\n")
}

# Create the matrix with NA values
frame.mat <- matrix(data = NA, nrow = Nrow, ncol = length(samples))
rownames(frame.mat) <- seq(1, Nrow, 1)
colnames(frame.mat) <- samples

# Fill in the matrix with data from FrameCount_List_df
for (s in samples) {
  rowVal <- as.numeric(rownames(FrameCount_List_df[[s]]))
  frame.mat[rowVal, s] <- FrameCount_List_df[[s]]$f.change
}


group_condition <- lapply(strsplit(samples, "_"), function(x) (x[1:2])) #it depends on the sample name
group_condition <- sapply(group_condition, paste, collapse = "_") 


diff_frame <- function(control_ind, sample_ind,...){
  c.df <- na.omit(as.data.frame(frame.mat[,control_ind]))
  a.df <- na.omit(as.data.frame(frame.mat[,sample_ind]))
  # check if there are replicates or not
  if (length(control_ind) == 1 & length(sample_ind) == 1) {
    diff <- merge(as.data.frame(c.df), as.data.frame(a.df), by = 0)
  } else {
    diff <- merge(as.data.frame(rowMeans(c.df)), as.data.frame(rowMeans(a.df)), by = 0)
  }
  diff <- column_to_rownames(diff, var = "Row.names")
  colnames(diff) <- c("F.F1_aver","F.F0_aver")
  diff <- diff %>% mutate("diff" = F.F1_aver - F.F0_aver) %>% na.omit()
  return(diff)
}

rank_diff <- function(df){
  diff_select <- df %>% arrange(desc(diff)) 
  diff.ind <- rownames(diff_select)
  return(diff_select)
}

TwoGroup_Compare <- function(control,sample,...){
  c.ind <- which(grepl(control, samples))
  a.ind <- which(grepl(sample, samples)) 
  diff<- diff_frame(c.ind,a.ind)
  return(diff)
}

control = "YPD_ctrl"
sample = "CSM"
Frameshift_df <- TwoGroup_Compare(control,sample)
#saveRDS(Frameshift_df,file.path(rds,paste(whichsample,"Frameshift_df.rds",sep = "_")))
