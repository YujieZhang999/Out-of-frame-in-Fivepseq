# # -----------------------------------------------------------------------
library("EnhancedVolcano")
library(tidyverse)
library(hrbrthemes)
library(viridis)

# # -----------------------------------------------------------------------
# Analyse the normal frame MS data to see the change at protein level when change to CSM
# # -----------------------------------------------------------------------
MS_data <- read_delim(file.path(input,"MS_canonical_protein.txt"), 
                      +     delim = "\t", escape_double = FALSE, 
                      +     trim_ws = TRUE)
MS_data <- as.data.frame(MS_data)
rownames(MS_data) <- MS_data$Genes


# # -----------------------------------------------------------------------
# Check protein abuandance of high and low frameshifted genes
# # -----------------------------------------------------------------------
whichsample = "CSM_ctrl"
Frameshift_df <- readRDS(file.path(rds,paste(whichsample,"Frameshift_df.rds",sep = "_")))

Frameshift_df$gene <- rownames(Frameshift_df)
Degradation_DIFF_high <- Frameshift_df %>% 
                        filter(F.F1_aver >= 0.2 & F.F0_aver <= -0.2) %>% 
                        mutate(group = "high")) 

#
Degradation_DIFF_low <- Frameshift_df %>% 
                        filter((F.F1_aver > 0 & F.F0_aver > 0 )) %>% 
                        mutate(group = "low")) 

#
Degradation_DIFF_all <- Frameshift_df %>% mutate(group = "all")) 

# # change the gene name
MS_data_gene = bitr(MS_data$UniProtIds, fromType="UNIPROT", toType=c("ORF"), OrgDb="org.Sc.sgd.db")
MS_data_2 <- MS_data %>% left_join(MS_data_gene, by =c("UniProtIds"="UNIPROT"))

Degradation_DIFF_high_pro <- Degradation_DIFF_high %>% 
                             left_join(MS_data_2, by = c("gene"="ORF"))  %>% 
                             filter(Qvalue <0.01) %>% 
                             na.omit()


Degradation_DIFF_low_pro <- Degradation_DIFF_low %>% 
                            left_join(MS_data_2, by = c("gene"="ORF")) %>% 
                            filter(Qvalue <0.01) %>% 
                            na.omit()


Degradation_DIFF_all_pro <- Degradation_DIFF_all %>% 
                            left_join(MS_data_2, by = c("gene"="ORF")) %>% 
                            filter(Qvalue <0.01) %>% 
                            na.omit()


Degradation_DIFF_protein <- rbind(Degradation_DIFF_high_pro,Degradation_DIFF_low_pro,Degradation_DIFF_all_pro)
#saveRDS(Degradation_DIFF_protein, file.path(rds,"Fig5_Degradation_DIFF_protein.rds"))






