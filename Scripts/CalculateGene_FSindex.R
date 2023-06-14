
# # -----------------------------------------------------------------------
# # Input data is FrameCount_List_df from ReadFivepseq_Frame_individual.R script
# # Extract f.change from list
# # for each sample, save the f.change in the matrix list
# # -----------------------------------------------------------------------
FrameCount_List_df <- readRDS(file.path(rds,paste(whichsample,"IndiV_FrameCount_list.rds",sep = "_")))
samples = names(FrameCount_List_df)

# Read transcript_descriptors.txt and check the dimentionality
t.des.name <- "transcript_descriptors.txt"
t.des.files <- list.files(fivepseq_count, t.des.name, full.names = T, recursive = T)
t.des <- read.table(t.des.files[[1]], sep = "\t", header = T)
Nrow <- dim(t.des)[1]

# Define expected values from transcript_descriptor file
Test_dim <- perform_checkpoint_test(Nrow, species_name="yeast")
perform_checkpoint_test <- function(Nrow, species_name) {
  # Define expected number of rows for each organism
  expected_rows <- list(
    yeast = 6600,
    human = 82724,
    bsubs = 4325,
    lpla = 3063
  )
  # Get the expected number of rows for the specified species
  expected_rows_species <- expected_rows[[species_name]]
  # Checkpoint test
  if (Nrow == expected_rows_species) {
    cat(paste0("Passed ", species_name, " checkpoint test!\n"))
  } else {
    cat(paste0("Failed ", species_name, " checkpoint test. Nrow value does not match expected value.\n"))
  }
}

# # -----------------------------------------------------------------------
# # Compare the frameshift index between rich and poor media
# # -----------------------------------------------------------------------
frame.mat <- matrix(data = NA, nrow = Nrow, ncol = length(samples))
rownames(frame.mat) <- seq(1, Nrow, 1)
colnames(frame.mat) <- samples

for (s in samples) {
  rowVal <- as.numeric(rownames(FrameCount_List_df[[s]]))
  frame.mat[rowVal, s] <- FrameCount_List_df[[s]]$f.change
}

group_condition <- lapply(strsplit(samples, "_"), function(x) (x[1:2])) #it depends on the sample name
group_condition <- sapply(group_condition, paste, collapse = "_") 

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
  diff <- diff %>% mutate("diff" = F.fs_rich - F.fs_poor) %>% na.omit()
  return(diff)
}

TwoGroup_Compare <- function(control,sample,FS.df){
  c.ind <- which(grepl(control, samples))
  a.ind <- which(grepl(sample, samples)) 
  diff <- diff_frame(c.ind,a.ind,FS.df)
  return(diff)
}

control <- "YPD_ctrl"
sample <- "CSM_ctrl"
FS.df <- frame.mat
Frameshift_df <- TwoGroup_Compare(control,sample,FS.df)
#saveRDS(Frameshift_df,file.path(rds,paste(whichsample,"Frameshift_df.rds",sep = "_")))
