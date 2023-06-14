
library(Biostrings)
library(seqinr)
library(coRdon)
library(dplyr)

# # -----------------------------------------------------------------------
# # CALCULATE CODON COMPOSITION OF THE FRAMESHIFT GENES IN FIVEPSEQ DATA
# # -----------------------------------------------------------------------
calculate_codon_ratio <- function(rds, frameshift_index, transcript_assembly_path) {
  transcript_assembly_path <- read.delim(transcript_assembly_path)
  transcript_sequences <- DNAStringSet(transcript_assembly_path$type)
  transcript_sequences$name <- transcript_sequences$ID
  transcript_codon <- codonTable(transcript_sequences)
  transcript_codon_count <- as.data.frame(transcript_codon@counts)
  transcript_codon_count$gene <- seq(0,length(transcript_sequences)-1,by =1) 
  row_sums <- apply(transcript_codon_count[,-65], 1, sum)
  codon_ratio <- (transcript_codon_count[,-65] / row_sums) * 100
  codon_ratio$gene <- as.character(transcript_codon_count$gene)
  # Load the data frame and join with the codon ratio data frame
  frameshift_index <- readRDS(frameshift_index)
  frameshift_index$gene <- rownames(frameshift_index)
  frameshift_index <- frameshift_index %>% na.omit()
  frameshift_codon <- codon_ratio %>% 
                      dplyr::left_join(frameshift_index, by = "gene") %>% 
                      na.omit() %>% 
                      dplyr::filter((fs_rich > 0) & (fs_poor < 0)) # select frameshifted genes 
  return(frameshift_codon)
}

# Input the path of the working directory
whichsample <- "human_riboseq"
whichspecies <- "h.sapiens"
whichcondition <- "293t_rpC_leu_3h"
rds <- file.path(workdir,"rds")
path <- file.path(workdir,paste(whichsample,"fivepseq",sep = "_")) # The path of the fivepseq output
transcript_assembly_path <- file.path(path, "transcript_assembly.txt")
frameshift_index <- file.path(rds,paste(whichcondition,"Frameshift_df.rds",sep = "_")) # The path of the frameshift index from CalculateGene_FSindex.R

FS_codon <- calculate_codon_ratio(rds, frameshift_index, transcript_assembly_path)

# # -----------------------------------------------------------------------
# # SET THE PERMUTATION TO CALCULATE THE CONFIDENTIAL INTERVAL
# # -----------------------------------------------------------------------
permutation_100 <- function(codon_df,sample_name,conf_level){
  n_permutations <- 100
  # Create an empty list to store the permuted data frames
  permuted_dfs <- vector("list", n_permutations)
  # Loop through the number of permutations and permute the 'x' column
  for (i in 1:n_permutations) {
    permuted_col <- sample(codon_df$diff) # Permute the column 'x'
    permuted_df <- codon_df # Create a copy of the original data frame
    permuted_df$x <- permuted_col # Replace the 'x' column with the permuted values
    permuted_dfs[[i]] <- permuted_df # Save the permuted data frame in the list
  }

  # Create an empty list to store the output matrices
  output_list <- vector("list", length(permuted_dfs))
  # Loop over the list of data frames
  for (i in seq_along(permuted_dfs)) {
    codon_cor.mat <- matrix(data = 0, nrow = 64, ncol = 1)
    rownames(codon_cor.mat) <- colnames(permuted_dfs[[i]])[1:64]
    colnames(codon_cor.mat)[1] <- "cor"
    # Loop over the codons and calculate the correlations
    for (j in seq(1:64)){
      dat <- permuted_dfs[[i]][,c(j,65,69)]
      codon <- colnames(dat)[1]
      corr <- cor(dat[,3], dat[,1], method="pearson")
      codon_cor.mat[codon, 1] <- corr
    }
    # Convert the output matrix to a data frame and store it in the output list
    output_list[[i]] <- as.data.frame(codon_cor.mat)
  }
  output_df <- do.call(cbind,output_list)
  # Define a function to calculate the confidence interval for a row
  calc_ci <- function(row) {
    n <- length(row)
    mean_val <- mean(row)
    std_err <- sd(row)/sqrt(n)
    q <- qnorm((1+conf_level)/2)
    lower_ci <- mean_val - q*std_err
    upper_ci <- mean_val + q*std_err # nolint
    return(c(lower_ci, upper_ci))
  }
  # Apply the function to each row of the data frame
  ci_df <- apply(output_df, 1, calc_ci)
  ci_df <- t(ci_df)
  colnames(ci_df) <- c("ci_min", "ci_max")
  saveRDS(ci_df, file.path(rds,paste(whichsample,"codon_permution.rds",sep = "_")))
  return(ci_df)
}

conf_level <- 0.95
confidential_interval_df <- permutation_100(FS_codon, whichsample, conf_level)

# # -----------------------------------------------------------------------
# # CALCULATE THE CORRELATION BETWEEN THE FRAMESHIFT INDEX AND THE CODON RATIO
# # -----------------------------------------------------------------------
codon_fs_correlation <- function(codon_df,ci_df){
  codon_cor.mat <- matrix(data = 0, nrow = 64, ncol = 1)
  rownames(codon_cor.mat) <- colnames(codon_df)[1:64]
  colnames(codon_cor.mat)[1] <- "cor"

  for (i in seq(1:64)){
    dat <- codon_df[,c(i,65,68)] # correlated with the frameshift index (diff)
    codon <- colnames(dat)[1]
    corr <- cor(dat[,3], dat[,1], method="pearson")
    codon_cor.mat[codon, 1] <- corr
  }
  codon_cor.mat <- as.data.frame(codon_cor.mat)
  codon_cor.mat_merge <- cbind(codon_cor.mat,ci_df)
  codon_cor.mat_merge$codons <- rownames(as.data.frame(codon_cor.mat_merge))
  return(codon_cor.mat_merge)
}

codon_cor <- codon_fs_correlation(FS_codon, confidential_interval_df)

# # -----------------------------------------------------------------------
# # MERGE WITH THE CODON OPTIMALITY or tRNA ADAPTATION INDEX FROM DIFFERENT SPECIES
# # -----------------------------------------------------------------------
merge_plot_codon_frequency <- function(CodonIndex, codon_df,whichsample) {
  CodoIndex_df <- read_excel(CodonIndex)
  codon_cor.mat <- codon_df %>% dplyr::left_join(CodoIndex_df,by = "codons") %>% na.omit()
  codon_cor.mat$AA <- paste(codon_cor.mat$AA,codon_cor.mat$codons,sep = "-")
  codon_cor.mat_reorder <- codon_cor.mat %>% dplyr::arrange(desc(cor)) 
  codon_cor.mat_reorder$AA <- factor(codon_cor.mat_reorder$AA, levels = codon_cor.mat_reorder$AA)
  saveRDS(codon_cor.mat_reorder,file.path(rds,paste(whichsample, "_codon_frequency_barplot.rds", sep="_"))
} # Use the saved data for plot Fig6

CodonIndex <- file.path(path, paste(whichspecies,"tAI.xlsx",sep = "_")) # Download for different species
merge_plot_codon_frequency(CodonIndex,codon_cor,whichsample)