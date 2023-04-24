library("LSD")
source(file.path(source_code,"LSD.heatscatter.R"))
# # -----------------------------------------------------------------------
# This is to plot density plot using LSD package
DIFF <- readRDS(file.path(rds,paste(whichsample,"YPD_CSM_ctrl_merge_diff_mean.rds",sep = "_")))

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



















