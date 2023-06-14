
library(ggthemes)
library("tidyverse")
library(dplyr)
library(ggrepel) 
options(stringsAsFactors = FALSE)
library(minpack.lm)

slamdir<- "out-of-frame/slamseq"
conversiondir <- file.path(slamdir,"slam_seq_dataset/slamseq_tc1")
rds <- file.path(conversiondir,"rds")
plotdir <- file.path(slamdir,"plot")

# # -----------------------------------------------------------------------
# Boxplot for the comparation for different group of genes for individual samples
# # -----------------------------------------------------------------------
Degradation_df <- readRDS(file.path(rds,"Degradation_df_coding.rds")) # no CUTS and XUTS included
Degradation_df$condition <- sapply(strsplit(Degradation_df$group, "_"), 
                                              function(x) {o <- paste(x[c(1,2)], collapse = "_")})

level_order <- factor(Degradation_df$group, level = c('by_ypd_r1','by_ypd_r2','by_ypd_r3',
                                                  'upf1_ypd_r1','upf1_ypd_r2','upf1_ypd_r3',
                                                  'by_csm_r1', 'by_csm_r2','by_csm_r3',
                                                  'upf1_csm_r1','upf1_csm_r2','upf1_csm_r3')


p <- ggplot(Degradation_df, aes(x = level_order, y = DR, color = condition)) +
  ggforce::geom_sina(alpha = 0.1) + 
  ylim(c(0,7.5)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  stat_summary(fun = median, color = "black") 
ggsave(filename = file.path(plotdir,"Degradation_boxplot_individuals.pdf"), width = 30, height = 30, units = "cm")


# # -----------------------------------------------------------------------
# Boxplot for the comparation for different group of genes for merged replicates 
# # -----------------------------------------------------------------------
Degradation_df_coding <- readRDS(file.path(rds,"Degradation_df_coding.rds")) # save from SLAMseq_halflife.R script
level_order <- factor(Degradation_df_coding$group, level = c('by_ypd','upf1_ypd','by_csm','upf1_csm'))


# Use the code above for the boxplot
p <- ggplot(Degradation_df_coding, aes(x = level_order, y = mean_all, color = group)) +
  ggforce::geom_sina(alpha = 0.1) + 
  ylim(c(0,7.5)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  stat_summary(fun = median, color = "black") 
ggsave(filename = file.path(plotdir,"Degradation_boxplot_merge.pdf"), width = 30, height = 30, units = "cm")


# # -----------------------------------------------------------------------
# GO enrichment for NMD sensitive genes
# # -----------------------------------------------------------------------
Degradation_merge <- readRDS(file.path(rds,"Degradation_merge.rds"))
NMD_sensitive <- readRDS(file.path(rds,"NMD_sensitive.rds"))


performGOEnrichment <- function(gene, gene_merge, ont) {
  eg <- bitr(as.factor(gene$gene), fromType = "ORF", toType = c("ENTREZID"), OrgDb = "org.Sc.sgd.db")
  eg_r <- bitr(as.factor(gene_merge$gene), fromType = "ORF", toType = c("ENTREZID"), OrgDb = "org.Sc.sgd.db")
  
  ego2 <- enrichGO(gene = eg$ENTREZID,
                   OrgDb = org.Sc.sgd.db,
                   keyType = 'ENTREZID',
                   universe = eg_r$ENTREZID,
                   ont = ont,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
  
  egobp2 <- simplify(ego2, cutoff = 0.6, by = "p.adjust", select_fun = min)
  
  df_ego2 <- data.frame(ego2)
  df_ego2_filter <- df_ego2[df_ego2[6] <= 0.05, ][, c(2, 6, 9)]
  df_ego2_filter$GO <- rownames(df_ego2_filter)
  df_ego2_filter$description <- paste(df_ego2_filter$Description, "(", df_ego2_filter$GO, ")" )
  rownames(df_ego2_filter) <- df_ego2_filter$description
  df_ego2_filter <- df_ego2_filter[1:5,2:3]
  return(df_ego2_filter)
}

# Perform enrichment for BP (Biological Process)
df_ego2_filter_cc <- performGOEnrichment(NMD_sensitive, Degradation_merge, "BP")
# Perform enrichment for MF (Molecular Function)
df_ego2_filter_mf <- performGOEnrichment(NMD_sensitive, Degradation_merge, "MF")
# Perform enrichment for CC (Cellular Component)
df_ego2_filter_cc <- performGOEnrichment(NMD_sensitive, Degradation_merge, "CC")

NMD_sensitive_GO <- rbind(df_ego2_filter_cc,df_ego2_filter_mf,df_ego2_filter_cc)
#write.csv(NMD_sensitive_GO, file.path(output,"NMD_sensitive_GO.csv"), row.names=TRUE,quote = FALSE)

NMD_sensitive_GO %>%
  ggplot(aes(x=rownames(NMD_sensitive), y=log10(p.adjust), fill=group)) +
  geom_col(position="dodge") + 
  coord_flip() + 
  labs(x="NMD_sensitive", y="log10(p.adjust)")
ggsave(filename = file.path(plotdir,"NMD_sensitive_GO.pdf"), width = 20, height = 20, units = "cm")




