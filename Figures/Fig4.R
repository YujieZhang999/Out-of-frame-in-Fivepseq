library(dplyr) 
library(ggplot2)
library(ggpubr)
library(pheatmap)
# # -----------------------------------------------------------------------
# Heatmap for codon frameshift index, related to Figure 4 and Figure S6
# # -----------------------------------------------------------------------
# Read output from CalculateCodon_FSindex.R
Frameshift_codon <- readRDS(file.path(rds,"Frameshift_codon.rds"))
whichsample <- "CSM_ctrl"
breaksList = seq(-2, 2, by = 0.2)
pdf(file.path(plotdir, paste(whichsample,"_heatmap_scale_codon.pdf",sep="")),width = 3.8, height = 9)
pheatmap(data.matrix(riboseq[,c(3,1)]),
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         #color=colorRampPalette(c("navy", "white", "red"))(20),
         color=colorRampPalette(c("navy", "white", "red"))(length(breaksList)),
         annotation_row = riboseq[,c(4,5,6)],
         breaks = breaksList, # set scales
         show_rownames = TRUE,fontsize_row = 7)
dev.off()


# # -----------------------------------------------------------------------
# plot codon correlation between fivepseq and riboseq
# # -----------------------------------------------------------------------
Riboseq_codon <- readRDS(file.path(rds,"Riboseq_codon.rds"))

CSM_ctrl <- CSM_ctrl %>% rename(diff_fivepseq = diff) %>% select(diff_fivepseq,AA,Coller,optimal)
Riboseq_codon <- Riboseq_codon %>% rename(diff_riboseq = diff) %>% select(diff_riboseq,AA)
merge_codon <- CSM_ctrl %>% 
                left_join(riboseq, by = "AA") %>% 


ggp <- ggplot(merge_codon, aes(diff_fivepseq, diff_riboseq,label = AA)) +    
  geom_point() +
  theme_bw() +
  geom_smooth(method = lm, fill = "lightgray") +
  stat_cor(label.x = 0.45, label.y = 1.0, method = "pearson",size= 4)+
  geom_text_repel(aes(size = 3, color=factor(group)))
ggp
ggsave(filename = file.path(plotdir,paste(whichsample,"_codon_correlated_fivepseq",".pdf")), width = 20, height = 20, units = "cm")


# # -----------------------------------------------------------------------
# plot AA addition with pairs
# # -----------------------------------------------------------------------
whichsample = 'CSM_AA_new'
CSM_AA_new <- readRDS(file.path(rds,"CSM_AA_new_codon.rds"))
CSM_AA_new$codon <- rownames(CSM_AA_new)

CSM_AA_new_long <- gather(CSM_AA_new, group, FS_index, A:S, factor_key=TRUE)
CSM_AA_new_A <-CSM_AA_new_long %>% dplyr::filter(group=="S" | group=="C") 
#CSM_AA_new_A$target_codon <- CSM_AA_new_A$codon %in% c("SER_TCT","SER_TCC","SER_TCA","SER_TCG","SER_AGT","SER_AGC")
CSM_AA_new_A$target_codon <- CSM_AA_new_A$codon %in% c("PRO_CCT","PRO_CCC","PRO_CCA","PRO_CCG")


CSM_AA_new_A$group <- factor(CSM_AA_new_A$group, levels = c("C", "S"))
p <- ggpaired(CSM_AA_new_A, x = "group", y = "FS_index",
              color = "group", palette = "jco", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = FALSE)
p + stat_compare_means(label = "p.format", paired = TRUE)
ggsave(filename = file.path(plotdir,paste(whichsample,"_codon_connected_Serine",".pdf")), width = 20, height = 20, units = "cm")


ggplot(CSM_AA_new_A,aes(x = group, y = FS_index)) +
  geom_boxplot(aes(fill = group), alpha = 0.2, col = "grey") +
  geom_point() +
  geom_line(aes(group = codon, col = target_codon)) +
  scale_colour_manual(values = c("black", "red"))
ggsave(filename = file.path(plotdir,paste(whichsample,"_codon_connected_Serine_color",".pdf")), width = 20, height = 20, units = "cm")





