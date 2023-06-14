library("EnhancedVolcano")
library(tidyverse)
library(hrbrthemes)
library(viridis)
# # -----------------------------------------------------------------------
# # Analyse the normal frame MS data to see the change at protein level when change to CSM
# # -----------------------------------------------------------------------
MS_data <- read_delim(file.path(input,"MS_canonical_protein.txt"), 
                      +     delim = "\t", escape_double = FALSE, 
                      +     trim_ws = TRUE)
MS_data <- as.data.frame(MS_data)
EnhancedVolcano(MS_data,
                lab = rownames(MS_data),
                x = 'AVG Log2 Ratio',
                y = 'Qvalue',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)
ggsave(file.path(plotdir,"Fig5_volcano_MS.pdf"), width = 30, height = 30, units = "cm")

# # -----------------------------------------------------------------------
# Analyse the normal frame MS data to see the change at protein level when change to CSM
# # -----------------------------------------------------------------------
Degradation_DIFF_protein <- readRDS(Degradation_DIFF_protein, file.path(rds,"Fig5_Degradation_DIFF_protein.rds"))
ggplot(data = Degradation_DIFF_protein,aes(x = group, y = `AVG Log2 Ratio`, fill = group))+
  scale_fill_viridis_d( option = "D")+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="gray30",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="gray30",alpha=0.5)+ 
  stat_summary(fun = median, width = 1.15, linetype = "solid")+
  theme_pubr()
ggsave(file.path(plotdir,"Fig5_protein_high_low_abundance.pdf"), width = 15, height = 30, units = "cm")

# # -----------------------------------------------------------------------
# Protein enrichment GSEA
# # -----------------------------------------------------------------------
GSEA_enrichment.Component <- read_delim(file.path(input,"GSEA_enrichment.Component.tsv")
GSEA_enrichment.Component_2 <- GSEA_enrichment.Component %>% filter(direction=="bottom") %>% arrange(desc(false.discovery.rate))

GSEA_enrichment.Component_2$term.description <- factor(GSEA_enrichment.Component_2$term.description, levels=unique(GSEA_enrichment.Component_2$term.description))
ggplot(GSEA_enrichment.Component_2, aes(x = `term.description`, y = log10(`false.discovery.rate`))) +
  geom_bar(stat = "identity",
           show.legend = TRUE,
           aes(fill = `enrichment.score`),  # Background color
           color = "gray30") + # Border color
  facet_grid(. ~ direction,scales="free") +
  xlab("Group") +
  ylab("Value") +
  coord_flip() +
  theme_minimal() +
  theme( axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  scale_fill_gradient2()
ggsave(filename = file.path(plotdir,"Fig5_Protein_GSEA_comp.pdf"), width = 30, height = 30, units = "cm")




