
# # -----------------------------------------------------------------------
# This is to calculate the frame preference for each individual samples without
# merge the index together
# # -----------------------------------------------------------------------

#workdir <- "/Users/yujiezhang/OneDrive - KI.SE/2021_slamseq_projects/2021_slamseq_for_FS"
#fivepseq_count <- file.path(workdir,"fivepseq_CSM_AA")
#rds <- file.path(workdir,"rds_repeat")
#plotdir <- file.path("/Users/yujiezhang/OneDrive - KI.SE/2021_slamseq_projects/2021_slamseq_for_FS/plot","Stress")
#dataDir <- file.path("/Users/yujiezhang/OneDrive - KI.SE/2021_slamseq_projects/2021_slamseq_for_FS","stress")
#datadir <- file.path(workdir,"datadir")
#whichsample ="CSM_AA"
# Load libraries

library("tidyverse")
library(dplyr)
library("zoo") # Using index
library("cowplot")
library("coRdon")
library("RGenetics")
library("ggpubr")
library("ggrepel")
opts <- options(stringsAsFactors = F)
## data input: generates a list with samples, and keeps the frame stats, transcript descriptors and libsize for each sample. 

data.name <- "data_summary.txt"
data.files <- list.files(fivepseq_count, data.name, full.names = T, recursive = T)
samples <- unlist(lapply(data.files, function(x) {
  basename(dirname(x))
}))

frame.name <- "transcript_frame_prefs.txt"
frame.files <- list.files(fivepseq_count, frame.name, full.names = T, recursive = T)
mydata <- list()
for (frame.f in frame.files) {
  frame <- read.table(frame.f, sep = "\t", header = T)
  name <- basename(dirname(frame.f))
  mydata[[name]][["frame"]] <- frame
}

t.des.name <- "transcript_descriptors.txt"
t.des.files <- list.files(fivepseq_count, t.des.name, full.names = T, recursive = T)
for (t.des.f in t.des.files) {
  t.des <- read.table(t.des.f, sep = "\t", header = T)
  name <- basename(dirname(t.des.f))
  mydata[[name]][["libsize"]] <- sum(t.des$NumOfReads)
  t.des <- t.des %>% mutate(rpm = 10^6 * t.des$NumOfReads / mydata[[name]]$libsize )
  mydata[[name]][["t.des"]] <- t.des
}

## count positions are extended 100 nt from the start and the end.
full_length.name <- "counts_FULL_LENGTH.txt"
count.length.files <- list.files(fivepseq_count, full_length.name, full.names = T, recursive = T)
for (count.f in count.length.files) {
  l <- readLines(count.f)
  name <- basename(dirname(count.f))
  counts <- strsplit(l, "\t") %>%  map(.,as.numeric)
  mydata[[name]][["count.f"]] <- counts
}

t.ass.name <- "transcript_assembly.txt"
t.ass.files <- list.files(fivepseq_count, t.ass.name, full.names = T, recursive = T)
for (t.ass.f in t.ass.files) {
  t.ass <- read.table(file.path(workdir, "transcript_assembly.txt"), sep = "\t", header = T)
  t.ass$ID <- sapply(strsplit(t.ass$ID, "\\:| |\\_| "), function(x) (x[2]))
  name <- basename(dirname(t.ass.f))
  t.ass <- t.ass %>% mutate(cds_start = cds_start + 1) # This is to match the 5pseq reference with database
  mydata[[name]][["t.ass"]] <- t.ass
}
#saveRDS(t.ass,file.path(workdir,"t.ass.rds"))
# filter low expressed genes: generated and index,which indicates high coverage genes
# if only read the output from last step
# retain.ind <- c()
# rpm_threshold <- 20
# count_threshold <- 80
# pos_threshold <- 40
# gene_length <- 400 ## select genes longer than 400-bp 
# # For low sequenced samples

retain.ind <- c()
rpm_threshold <- 10
count_threshold <- 30
pos_threshold <- 20
gene_length <- 400

for (s in samples) {
  t.des <- mydata[[s]][["t.des"]]
  pass.ind <- which(10^6 * t.des$NumOfReads / mydata[[s]]$libsize >= rpm_threshold &
                      t.des$NumOfReads >= count_threshold &
                      t.des$NumOfMapPositions >= pos_threshold &
                     # t.des$X3nt == "True" &  These parameters are for human, becasue they have too entires
                    #  t.des$start == "ATG"  &
                     # t.des$stop %in% c("TAG","TGA","TAA") &
                      t.des$len >= gene_length)
  
  retain.ind <- c(retain.ind, pass.ind)
  retain.ind <- intersect(pass.ind, retain.ind)
  retain.ind <- unique(retain.ind)
  mydata[[s]][["pass.ind"]] <- retain.ind
}


# Save a count without extension of start and end
sublist_without_Ext <- function(whichlist = x1) {
  myfilter <- whichlist[["count.f"]] %>% map_int(~length(.x))
  mydf <- map2(whichlist[["count.f"]],myfilter, ~.x[101:(.y-100)]) ## Extract from 101 to the end of gene (because 5pseq count 100 from the start and the end)
  newname_df <- paste0("ORF_", "count")
  whichlist[[newname_df]] <- mydf
  myframe <- map(mydf,~index(.x)%%3) ## Calculate frame for each gene
  newname_frame <- paste0("ORF_", "frame")
  whichlist[[newname_frame]] <- myframe
  return(whichlist)
}
mydata <- map(mydata, sublist_without_Ext)

## Extract ORF frame and count using the sub genes list
SubFrame_list <- map(mydata, ~.x$ORF_frame[.x$pass.ind])
SubCount_list <- map(mydata, ~.x$ORF_count[.x$pass.ind])

## Exclude the first 50bp and last 50bp (exclude non-frame-shift region)
SubList_Count_Exlude50 <- function(whichlist) {
  myfilter <- whichlist %>% map_int(~length(.x))
  mydf <- map2(whichlist, myfilter, ~.x[51:(.y-50)])
  return(mydf)
}

SubFrame_list_Exclude50 <- map(SubFrame_list, SubList_Count_Exlude50)
SubCount_list_Exclude50 <- map(SubCount_list, SubList_Count_Exlude50) 

SubCountFrame_list <- list("SubFrame_list"=SubFrame_list,
                           "SubCount_list"=SubCount_list,
                           "SubFrame_list_Exclude50"=SubFrame_list_Exclude50,
                           "SubCount_list_Exclude50"=SubCount_list_Exclude50)


## Calculate theboo gene specific frame coverage using the middle of the genes
CountFrame <- function(whichFrame = x1,whichCount = x2) {
  F0_index <- map(whichFrame, ~which(.x==0))
  F1_index <- map(whichFrame, ~which(.x==1))
  F2_index <- map(whichFrame, ~which(.x==2))
  F0_SumCount <- map2(whichCount,F0_index, ~sum(.x[.y]))
  F1_SumCount <- map2(whichCount,F1_index, ~sum(.x[.y]))
  F2_SumCount <- map2(whichCount,F2_index, ~sum(.x[.y]))
  
  FrameSum <- list(F0_TotolCount=unlist(F0_SumCount),
                   F1_TotolCount=unlist(F1_SumCount),
                   F2_TotolCount=unlist(F2_SumCount))
  return(FrameSum)
}

FrameCount_list <- map2(SubFrame_list_Exclude50,SubCount_list_Exclude50,CountFrame)
FrameCount_ORF_list <- map2(SubFrame_list,SubCount_list,CountFrame)#this is the control to check the frame 

# from list into data frame
FrameCount_df <- FrameCount_list %>% map(~ t(Reduce(rbind,.)))
FrameCount_ORF_df <- FrameCount_ORF_list %>% map(~ t(Reduce(rbind,.)))

# change colnames for all data frame
for (i in 1:length(FrameCount_df)){
  colnames(FrameCount_df[i][[1]]) <- c("F2","F0","F1")
  rownames(FrameCount_df[i][[1]]) <- mydata[[i]][["pass.ind"]]
}
for (i in 1:length(FrameCount_ORF_df)){
  colnames(FrameCount_ORF_df[i][[1]]) <- c("F2","F0","F1")
  rownames(FrameCount_ORF_df[i][[1]]) <- mydata[[i]][["pass.ind"]]
}

Calc_percent <- function(dataframe = x1){
  dataframe <- as.data.frame(dataframe)
  dataframe2 <- dataframe %>% dplyr :: mutate(total = (apply(dataframe, 1, sum)), 
                                              F0_per = F0 / total, F1_per = F1 / total, F2_per = F2 / total, 
                                              max_Frame = max.col(dataframe[,1:3]))
  dataframe2$Max_F[dataframe2$max_Frame == 3] <- 1
  dataframe2$Max_F[dataframe2$max_Frame == 2] <- 0
  dataframe2$Max_F[dataframe2$max_Frame == 1] <- 2
  dataframe3 <- dataframe2 %>% mutate(fpi=NA,
                                      fpi=ifelse(Max_F==1,log2(F1/((F2+F0)/2)),fpi),
                                      fpi=ifelse(Max_F==0,log2(F0/((F2+F1)/2)),fpi),
                                      fpi=ifelse(Max_F==2,log2(F2/((F1+F0)/2)),fpi), 
                                      f.change = log2(F1/F0)) %>% dplyr ::select(-max_Frame) # correct for riboseq is F0/F2
  return(dataframe3)
}

## Dataframe with frame counts, fpi and f.change (f.change = log2(F1/F0))
FrameCount_List_df <- map(FrameCount_df,Calc_percent) 
FrameCount_ORF_List_df <- map(FrameCount_ORF_df,Calc_percent) # This is the control

# This is to transform lists into a data frame for all samples (including fpi)
#FrameCount_df <- map_df(FrameCount_List_df, ~as.data.frame(.x), .id="Sample")
#FrameCount_ORF_df <- map_df(FrameCount_ORF_List_df, ~as.data.frame(.x), .id="Sample")


#saveRDS(FrameCount_df,file.path(rds,paste(whichsample,"_FrameCount_df.rds")))
saveRDS(FrameCount_List_df,file.path(rds,paste(whichsample,"IndiV_FrameCount_list.rds",sep = "_")))
saveRDS(FrameCount_ORF_List_df,file.path(rds,paste(whichsample,"IndiV_FrameCount_ORF_list.rds",sep = "_")))

