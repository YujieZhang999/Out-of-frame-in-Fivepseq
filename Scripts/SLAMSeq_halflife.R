
library("tidyverse")
library(dplyr)
library(ggrepel) 
options(stringsAsFactors = FALSE)
library(minpack.lm)

slamdir<- "out-of-frame/slamseq"
conversiondir <- file.path(slamdir,"slam_seq_dataset/slamseq_tc1")
rds <- file.path(conversiondir,"rds")
plotdir <- file.path(slamdir,"plot")

# Read tc counts reads from SLAM-DUNK
slam <- readRDS("NormTotal_OneConver.rds")
slamdata <- slam %>% map(select, gene_name, norma_tcCount) %>% imap(~ set_names(.x, c("gene_name", paste0("norma_tcCount_", .y))))

# # -----------------------------------------------------------------------
# Double check the list name if corresponding to the names in the rds data
# # -----------------------------------------------------------------------
slamlist <- list(
  by_csm_r1 = paste0("by_csm_r1_slam_t", c(0, 15, 30)),
  by_ypd_r1 = paste0("by_ypd_r1_slam_t", c(0, 15, 30)),
  upf1_csm_r1 = paste0("upf1_csm_r1_slam_t", c(0, 15, 30)),
  upf1_ypd_r1 = paste0("upf1_ypd_r1_slam_t", c(0, 15, 30)),
  by_csm_r2 = paste0("by_csm_r2_slam_t", c(0, 15, 30)),
  by_ypd_r2 = paste0("by_ypd_r2_slam_t", c(0, 15, 30)),
  upf1_csm_r2 = paste0("upf1_csm_r2_slam_t", c(0, 15, 30)),
  upf1_ypd_r2 = paste0("upf1_ypd_r2_slam_t", c(0, 15, 30)),
  by_csm_r3 = paste0("by_csm_r3_slam_t", c(0, 15, 30)),
  by_ypd_r3 = paste0("by_ypd_r3_slam_t", c(0, 15, 30)),
  upf1_csm_r3 = paste0("upf1_csm_r3_slam_t", c(0, 15, 30)),
  upf1_ypd_r3 = paste0("upf1_ypd_r3_slam_t", c(0, 15, 30))
)

slamnewdata <- map(slamlist,~ slamdata[.x]) %>%
               map(~ purrr::reduce(., left_join, by = "gene_name")) %>%
               map(~ filter_at(., vars(starts_with("norma")), all_vars(. > 20)))
glimpse(slamnewdata)

# # -----------------------------------------------------------------------
# Fit the Non-linear decay model equation from Herzog et al., 2017  
# # -----------------------------------------------------------------------
library(parallel)
tm <- c(0, 15, 30)
calc.decay.m1 <- function(y, tm) {
         # control <- list(maxiter=50, tol=1e-5)
          start <- list(#a=log(2)/10, # assume 10min half life
                        C=1,            
                        Plat = 0,
                        a=log(2)/10)
          fit <- nlsLM(y ~ Plat + (C-Plat)*exp(-a*(tm)), start=start, 
                       control = nls.lm.control(maxiter = 1000),
                       upper = c(1,0,Inf),
                       lower = c(1,0,0.001),
                       na.action = na.omit,trace = TRUE)
          C.hat <- coef(fit)[1]
          a.hat <- coef(fit)[3]
          Plat.hat <- coef(fit)[2]
          y.hat <- Plat.hat + (C.hat-Plat.hat) * exp(-a.hat * tm)
          resid <- sum((y-y.hat)**2)
          T.half.hat <- log(2)/coef(fit)["a"]
          list(a.hat=a.hat, C.hat=C.hat, y.hat=y.hat, T.half.hat=T.half.hat, fit=fit,
          resid=resid)
          }


HALF.LIFE <- list()
# calculate the half time for each strain/condition in the SAMPLE.SETS list
for(j in 1:length(slamnewdata)) {
  # extract the tpm matrix with the above column number
  X.analyze_sub <- slamnewdata[[j]]
  X.analyze <- subset(X.analyze_sub, (X.analyze_sub[,2] > X.analyze_sub[,3]) & (X.analyze_sub[,3] > X.analyze_sub[,4]))
  rn <- X.analyze %>% pull(gene_name)
  X.analyze_matrix <- X.analyze[,2:4]
  start <- proc.time()  # initial starting time for the whole process
  # using parallel computing for each gene/row in the matrix
  # nrow(X.analyze), is the number of the genes in total
  x <- mclapply(1:nrow(X.analyze_matrix), mc.cores=8, FUN=function(i) {
    #if(i %% 1000 == 0) print(
    #sprintf("[%d:%d] [%d:%d] Calculating half life %.2fs",j, length(SAMPLE.SETS), i, nrow(X.analyze_matrix), (proc.time()-start)[3]))
    Y.t <- as.numeric(X.analyze_matrix[i,])
    # Normalize to first time point.
    Y.t <- Y.t/Y.t[1]
    Y.t <- Y.t[Y.t >0 & !is.nan(Y.t)]
    m1 <- tryCatch({
    calc.decay.m1(Y.t, tm)
    }, error=function(err) {
    as.character(err)
    })
    x <- rep("", 5)
    if(is.character(m1)) {
    x[5] <- gsub("[\r\n]", "", m1)[1]
    } else {
    x <- c(m1$a.hat, m1$C.hat, m1$T.half.hat, m1$resid, "")
    }
    x
    })
  data.out <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  colnames(data.out) <- c("M1_a", "M1_C", "M1_T_half", "M1_residual", "M1 err")
  rownames(data.out) <- rn
  data.out <- data.out %>% as.data.frame %>% rownames_to_column(var = "gene_name") %>% tibble()
  HALF.LIFE[[names(slamnewdata)[j]]] <- data.out
}
halflife_coding <- HALF.LIFE %>% map(~ filter(., M1_T_half > 0)) %>% 
            map(filter, !grepl("CUT|SUT", gene_name))
  
#saveRDS(halflife_coding, file.path(rds,"halflife_OneConverNormTotal.rds")) # saved and for Fig3 plotting

# # -----------------------------------------------------------------------
# Calculate RNA degradation rate for individual samples
# # -----------------------------------------------------------------------
HL_slamseq <- readRDS(file.path(rds,"halflife_OneConverNormTotal.rds"))
samples <- names(HL_slamseq)

# Extracted relevant columns 
H_life <- map(HL_slamseq,~.x[ , c("gene_name", "M1_T_half")])

# Calculate degradation rate
DegRate <- map(H_life, ~ dplyr::mutate(.,M1_T_half=as.numeric(M1_T_half),DR = (0.693/(M1_T_half))*60))
DegRate <-  map(DegRate, ~.x[,c("gene_name","DR")]) 

DegRate <-  Map(function(x, n) setNames(x, c(names(x)[1], n)), DegRate, names(DegRate))
#save in Degradation_Rate.rds 

# # -----------------------------------------------------------------------
# Calculate RNA degradation rate by takeing the average for replicates
# # -----------------------------------------------------------------------
DegRate <- readRDS("Degradation_Rate.rds")
processData <- function(data, samples, group, output) {
  sample_list <- data[grep(samples, data)]
  result <- Reduce(function(x, y) merge(x, y, by = "gene_name"), sample_list)
  result <- result %>% mutate(mean_all = (DR.x + DR.y + DR) / 3)
  result$group <- group
  # write.xlsx(result, file.path(output, paste0(group, ".xlsx")), rowNames = TRUE)
  return(result)
}

by_csm <- processData(DegRate, "by_sc_", "by_csm", output)
upf1_csm <- processData(DegRate, "upf1_sc_", "upf1_csm", output)
by_ypd <- processData(DegRate, "by_ypd_", "by_ypd", output)
upf1_ypd <- processData(DegRate, "upf1_ypd_", "upf1_ypd", output)

degradation_merge <- rbind(by_csm,upf1_csm,by_ypd,upf1_ypd)
Degradation_df_coding <- degradation_merge %>% 
                        select(gene_name,mean_all,group) %>% 
                        map_df(Degradation_df_coding, ~as.data.frame(.x), .id="group")
#saveRDS(Degradation_df_coding,file.path(rds,"Degradation_df_coding.rds"))


# # -----------------------------------------------------------------------
# Degradation rate among all conditions
# # -----------------------------------------------------------------------
processDegradation <- function(Degradation_df, group1, group2) {
  Degradation_group1 <- Degradation_df %>% filter(group %in% group1)
  Degradation_group2 <- Degradation_df %>% filter(group %in% group2)
  colnames(Degradation_group1)[3] <- group1
  colnames(Degradation_group2)[3] <- group2

  Degradation_merge <- Degradation_group1 %>% left_join(Degradation_group2, by = "gene")
  Degradation_merge <- na.omit(Degradation_merge)
  Degradation_merge <- Degradation_merge %>% dplyr::select(gene, !!group1, !!group2) %>% dplyr::mutate(diff_df = !!group1 / !!group2)
  return(Degradation_merge)
}

Degradation_merge_csm <- processDegradation(Degradation_df, "by_csm", "upf1_csm")
Degradation_merge_ypd <- processDegradation(Degradation_df, "by_ypd", "upf1_ypd")

Degradation_merge <- Degradation_merge_csm %>% left_join(Degradation_merge_ypd, by ="gene") %>% na.omit() 
#saveRDS(Degradation_merge,file.path(rds,"Fig3_Degradation_merge.rds"))
#Degradation_merge <- readRDS(file.path(rds,"Degradation_merge.rds"))


# # -----------------------------------------------------------------------
# Define NMD sensitive genes
# Combine degradation rate and Frameshift index in YPD and CSM
# # -----------------------------------------------------------------------
NMD_sensitive <- Degradation_merge %>% filter(diff_df_ypd < 0.8 & diff_df_csm > 1.2)
#saveRDS(NMD_sensitive,file.path(rds,"NMD_sensitive.rds"))

Degradation_df$gene <- sapply(strsplit(Degradation_df$gene, ";"),function(x) {o <- x[1]})

whichsample = "CSM_ctrl"
Frameshift_df <- readRDS(file.path(rds,paste(whichsample,"Frameshift_df.rds",sep = "_")))
Frameshift_df$gene <- rownames(Frameshift_df)

DIFF_FS <- Frameshift_df %>%  filter(F.F1_aver > 0.2 & F.F0_aver < -0.2)

# To do the hypemetric analysis, I need to find number of successes k	, sample size s	, number of successes in the population M	, population size N	
NMD_sensitive_FS <- NMD_sensitive %>% left_join(DIFF_FS, by = "gene") %>% na.omit() # k = 246
NMD_sensitive_all <- NMD_sensitive %>% left_join(DIFF, by = "gene") %>% na.omit() # M = 712
NMD_resistance_all <- NMD_resistance %>% left_join(DIFF, by = "gene") %>% na.omit() # 545, N = 712 + 545 = 1257

total_FS <- DIFF_FS %>% left_join(Degradation_merge, by = "gene") %>% na.omit() %>% dim() # sample size s	= 433
dhyper(x = dim(NMD_sensitive_FS)[1], m = dim(NMD_sensitive_all)[1], n = dim(NMD_resistance_all)[1], k = total_FS[1],log = FALSE) # p = 0.04
#write.csv(DIFF_FS$gene, file.path(rds,"Frameshifted_gene.csv"), row.names=FALSE,quote = FALSE)