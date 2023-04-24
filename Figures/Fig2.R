
# Input data is FrameCount_List_df
# Extract f.change from list 
# for each sample, save the f.change in the matrix list

library("LSD")
source(file.path(source_code,"LSD.heatscatter.R"))

whichsample ="CSM_riboseq"
FrameCount_List_df <- readRDS(file.path(rds,paste(whichsample,"IndiV_FrameCount_list.rds",sep = "_")))
samples = names(FrameCount_List_df)


frame.mat <- matrix(data = NA, nrow = 6600, ncol = length(samples)) # yeast is 6600, human is 82724
rownames(frame.mat) <- seq(1,6600,1)
colnames(frame.mat) <- samples

for (s in samples) {
  rowVal <- as.numeric(rownames(FrameCount_List_df[[s]]))
  frame.mat[rowVal, s] <- FrameCount_List_df[[s]]$f.change
}
#saveRDS(frame.mat, file.path(rds,"human_riboseq_frame.mat.rds"))


group_condition <- lapply(strsplit(samples, "_"), function(x) (x[1:2])) #it depends on the sample name
group_condition <- sapply(group_condition, paste, collapse = "_") 


diff_frame <- function(control_ind, sample_ind,...){
  c.df <- as.data.frame(frame.mat[,control_ind])
  c.df <- na.omit(c.df)
  a.df <- as.data.frame(frame.mat[,sample_ind])
  a.df <- na.omit(a.df)
  
  # check if there are replicates or not
  if (length(control_ind) == 1 & length(sample_ind) == 1) {
    diff <- merge(as.data.frame(c.df), as.data.frame(a.df), by = 0)
  } else {
    diff <- merge(as.data.frame(rowMeans(c.df)), as.data.frame(rowMeans(a.df)), by = 0)
  }
  
  rownames(diff) <- NULL
  diff <- column_to_rownames(diff,var = "Row.names")
  colnames(diff) <- c("F.F1_aver","F.F0_aver")
  diff <- diff %>% mutate("diff" = F.F1_aver - F.F0_aver)
  diff <- na.omit(diff)
  return(diff)
}

rank_diff <- function(df){
  diff_select <- df %>% arrange(desc(diff)) 
  diff.ind <- rownames(diff_select)
  return(diff_select)
}

TwoGroup_Compare <- function(control,sample,...){
  c.ind <- which(grepl(control, samples))
  a.ind <- which(grepl(sample, samples)) # Which returns index
  diff_CSM <- diff_frame(c.ind,a.ind)
  rank_CSM <- rank_diff(diff_CSM)
  return(rank_CSM)
}

# # -----------------------------------------------------------------------
# This is to plot density plot using LSD package
# # -----------------------------------------------------------------------
# Specify the condition want to copare
DIFF <- TwoGroup_Compare("by4741_replete",unique(group_condition)[1])
#DIFF <- readRDS(file.path(rds,paste(whichsample,"YPD_CSM_ctrl_merge_diff_mean.rds",sep = "_")))
#saveRDS(DIFF,file.path(rds,paste(whichsample,"YPD_CSM_ctrl_merge_diff_median.rds",sep = "_")))

plot_heatscatter_density <- function(control,sample,df){
gg <- LSD.heatscatterpoints(df$F.F1_aver, df$F.F0_aver, xlab = control, ylab = sample, ggplot = TRUE)
gg + 
  geom_vline(xintercept = 0,linetype="dotted") +
  geom_hline(yintercept = 0,linetype="dotted") +
  
  xlim(c(-2,2)) + 
  ylim(c(-2,2)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

gg_dist_g1 = gghistogram(
  df, x = "F.F1_aver", y = "..density..",
  add = "median", rug = TRUE, 
  add_density = TRUE
)
  
gg_dist_g1 = gg_dist_g1 + ylab(control)

gg_dist_g2 = gghistogram(
  df, x = "F.F0_aver", y = "..density..",
  add = "median", rug = TRUE, 
  add_density = TRUE)
gg_dist_g2 = gg_dist_g2 + ylab(control)

# Flip axis of gg_dist_g2
gg_dist_g2 = gg_dist_g2 + coord_flip()
# Remove some duplicate axes
gg_dist_g1 = gg_dist_g1 + ylab(control) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) 

gg_dist_g2 = gg_dist_g2 + ylab(sample) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) 

# Modify margin c(top, right, bottom, left) to reduce the distance between plots
# and align G1 density with the scatterplot
gg_dist_g1 = gg_dist_g1 + theme(plot.margin = unit(c(0.5, 0, 0, 0.7), "cm"))
gg_scatter = gg + theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm")) +
              xlim(c(-2,2)) + 
              ylim(c(-2,2)) 
gg_dist_g2 = gg_dist_g2 + theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

# Combine all plots together and crush graph density with rel_heights
first_col = plot_grid(gg_dist_g1, gg_scatter, ncol = 1, rel_heights = c(1, 3))
second_col = plot_grid(NULL, gg_dist_g2, ncol = 1, rel_heights = c(1, 3))
perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1))
save_plot(filename = file.path(plotdir,paste(sample,control,"heatscatter_density.pdf",sep="_")), perfect, base_height=4)
}

for (i in 1:length(unique(group_condition))){
  df_Comp <- TwoGroup_Compare("YPD_ctrl",unique(group_condition)[i])
  plot_diff <- plot_heatscatter_density("YPD_ctrl",unique(group_condition)[i],df_Comp)
}


# # -------------Correlations with RNA featuures----------------------------------------------------------
# DIFF <- readRDS(file.path(rds,"CSM_ctrl_YPD_CSM_ctrl_merge_diff.rds"))
# 
# RNA_feature <- readRDS(file.path(rds,"RNA_feature.rds"))
# RNA_feature <- RNA_feature[c(1,3:5)]
# DIFF$Geneid <- rownames(DIFF)
# DIFF_features <- RNA_feature %>% map(~left_join(.x,DIFF,by="Geneid")) %>% map(~na.omit(.x))
#   
# 
# gg <- LSD.heatscatterpoints(DIFF_features[[4]]$diff, log2(DIFF_features[[4]]$median5), ggplot = TRUE)
# gg + 
#   xlim(c(-1,1)) + 
# #  ylim(c(0,5000)) +
#   theme_bw() + 
#   geom_smooth(method = lm, fill = "lightgray") +
#   stat_cor(label.x = -0.2, label.y = 0.25, method = "spearman",size=5) +
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"))
#  
# ggsave(filename = file.path(plotdir,paste(whichsample,"_correlation_diff_UTR5_mean",".pdf")), width = 20, height = 20, units = "cm")




















