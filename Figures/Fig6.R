library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library("tidyverse")
library(dplyr)
library(ggrepel)
options(stringsAsFactors = FALSE)

# # -----------------------------------------------------------------------
# # READ THE GENE FRAMESHIFT INDEX FROM CalculateGene_FSindex.R
# # -----------------------------------------------------------------------
whichsample <- "human_riboseq"
Frameshift_df <- readRDS(file.path(rds,paste(whichsample,"Frameshift_df.rds",sep = "_"))) 

# # -----------------------------------------------------------------------
# # PLOT FIG6, INCLUDING THE DENSITY PLOT AND CUMULATIVE PLOT
# # -----------------------------------------------------------------------
plot_distribution <- function(control, treatment, FS.df) {
  df_Comp <- TwoGroup_Compare(control, treatment, FS.df )

  data_long <- gather(df_Comp, sample, fs_index, fs_rich:fs_poor, factor_key=TRUE)
  phist <- gghistogram(
    data_long, x = "fs_index",
    add = "mean", rug = TRUE,
    fill = "sample", palette = c("#00AFBB", "#E7B800")
  )

  pdensity <- ggdensity(
    data_long, x = "fs_index",
    color= "sample", palette = c("#00AFBB", "#E7B800"),
    alpha = 0
  ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
    theme_half_open(11, rel_small = 1) +
    rremove("x.axis")+
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    rremove("legend")

  aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
  ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
  ggsave(filename = file.path(plotdir, paste(as.character(treatment),".pdf")), width = 20, height = 20, units = "cm")
  
  fs_test <- ks.test(df_Comp$fs_rich,df_Comp$fs_poor, exact=FALSE)
  ggplot(data_long, aes(x = fs_index)) +
    stat_ecdf(aes(color = sample,linetype = sample), geom = "step", size = 1.5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    annotate("text",x=0,y=1, label= paste("p-value =",fs_test$p.value)) +
    labs(y = "density") +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(filename = file.path(plotdir, paste(as.character(treatment),"ecdf.pdf")), width = 20, height = 20, units = "cm")
}

plot_distribution("293t_rpC_rich_3h","293t_rpC_leu_3h", Frameshift_df)
#plot_distribution("B.subtilis_rich","B.subtilis_low_nutrient")
#plot_distribution("L.plantarum_rich","L.plantarum_low_nutrient")
#plot_distribution("YPD","CSM")


# # -----------------------------------------------------------------------
# # PLOT CORRELATION OF FRAMESHIFT INDEX AND CODON FREQUENCY
# # -----------------------------------------------------------------------
whichsample <- "human_riboseq"
codon_cor.mat <- readRDS(file.path(rds,paste(whichsample, "_codon_frequency_barplot.rds", sep="_")))

t <- ggplot(codon_cor.mat, aes(x=AA, y=cor, fill=optimal)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=ci_min, ymax=ci_max), width=.2, position=position_dodge(.9)) +
  ylim(-0.2, 0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(plotdir, paste(whichsample, "codon_barplot_pearson.pdf",sep = "_")), width = 30, units = "cm")




