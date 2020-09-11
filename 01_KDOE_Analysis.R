require(tidyverse)
require(ggrepel)
require(reshape2)
require(viridis)
require(VennDiagram)

## FUNCTIONs
## Find pairs of duplicates (in the FullCO column, for example)
FindPairs = function (table, column) {
  Dups=which(duplicated(table[,column]))
  Pairs=list()
  for (i in rownames(table)[Dups]) {
    n=unlist(strsplit(i, split="-"))[1]
    lines=grep(n, rownames(table))
    Pairs[[n]]=table[lines,]
  }
  return(Pairs)
}  
## Check if all the duplicates found are only ALTA / ALTD events
CheckAlt = function (listofpairs) {
  n=names(listofpairs)
  if (is.null(n)) {return("EMPTY")
  } else { test=sapply(n, function(x) grepl("ALT",x))
  if (sum(test)==length(test)) {return("OK")
  } else {return("NO!")}
  }
}
## Count events (Wrap_VTS_Events function)
source("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/useful-R-scripts/VTS_Count_EventTypes.R")
## Separate events in a list
SplitList_inclMIC = function (input_list) {
  
  out <- list()
  CEx_inclMIC <- c("S","C1","C2","C3","MIC","ANN") 
  RIs <- c("IR-C","IR-S","IR")
  
  for(n in names(input_list)) {
    out[[n]] <- list(
      CEx = subset(input_list[[n]], COMPLEX %in% CEx_inclMIC),
      RI = subset(input_list[[n]], COMPLEX %in% RIs),
      Alt3 = subset(input_list[[n]], COMPLEX == "Alt3"),
      Alt5 = subset(input_list[[n]], COMPLEX == "Alt5")
    )}  
  return(out) }
## Plot scatterplots ddPSI (Plot_ScatterPlots_dPSI function)
source("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/R_Scripts/useful-R-scripts/Plot_Scatterplot_dPSI_values.R")





## IMPORT INCLUSION TABLES WITH ALL DiffAS
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/")
INCL <- list(
  KD=read.table("KD/DiffAS_dPSI10/INCLUSION_LEVELS_FULL-Mm218_withdPSI_9DiffAS.txt",header=T,sep="\t",fill=NA),
  OE=read.table("OE/DiffAS_dPSI10/INCLUSION_LEVELS_FULL-Mm218_withdPSI_3DiffAS.txt",header=T,sep="\t",fill=NA) 
)
CTRLs <- list(
  d12.d00_shSCR_all=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shSCR-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_Empty_all=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_Empty-with_dPSI.tab",header=T,sep="\t")
)
INCL <- lapply(INCL, function(x) {rownames(x)=x$EVENT; x})  
CTRLs <- lapply(CTRLs, function(x) {rownames(x)=x$EVENT; x})  

## IMPORT DIFFAS DATA FOR KD
VTS_KD_dPSI0 <- list(
  d12.d00_shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shSCR-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shC1=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shC5=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shU3=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shU4=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t"),
  d12_shC1.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
  d12_shC5.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
  d12_shU3.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
  d12_shU4.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t")
)
VTS_KD <- list(
  d12.d00_shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shSCR-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shC1=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shC5=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shU3=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shU4=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t"),
  d12_shC1.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
  d12_shC5.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
  d12_shU3.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
  d12_shU4.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t")
)
VTS_KD_dPSI0 <- lapply(VTS_KD_dPSI0, function(x) {rownames(x)=x$EVENT; x})  
VTS_KD <- lapply(VTS_KD, function(x) {rownames(x)=x$EVENT; x})  
lapply(VTS_KD_dPSI0, function(x) dim(x))  
lapply(VTS_KD, function(x) dim(x))  

## IMPORT DIFFAS DATA FOR OE
VTS_OE_dPSI0 <- list(
  d12.d00_Empty=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_TIA1=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
  d12_TIA1.Empty=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day12_Empty-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t")
)
VTS_OE <- list(
  d12.d00_Empty=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_Empty-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_TIA1=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
  d12_TIA1.Empty=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day12_Empty-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t")
)
VTS_OE_dPSI0 <- lapply(VTS_OE_dPSI0, function(x) {rownames(x)=x$EVENT; x})  
VTS_OE <- lapply(VTS_OE, function(x) {rownames(x)=x$EVENT; x})  


# Find duplicate lines with FullCO  (Alt3 and Alt5)
Pairs=lapply(VTS_KD, function(x) FindPairs(table=x, column = "FullCO"))  
summary(Pairs) #16 SPF45, 11 SR140, 14 CHERP, only 1 in UP_SPF45 and UP_SR140
# Check if all duplicates are ALT3 or ALT3
lapply(Pairs, function(x) CheckAlt(x))

# Find duplicate lines with FullCO  (Alt3 and Alt5)
Pairs=lapply(VTS_OE, function(x) FindPairs(table=x, column = "FullCO"))  
summary(Pairs) #16 SPF45, 11 SR140, 14 CHERP, only 1 in UP_SPF45 and UP_SR140
# Check if all duplicates are ALT3 or ALT3
lapply(Pairs, function(x) CheckAlt(x))

# Calculate how many events in each table
Numbers_KD=lapply(VTS_KD, function(x) summary(x[,"COMPLEX"])); Numbers_KD
Numbers_OE=lapply(VTS_OE, function(x) summary(x[,"COMPLEX"])); Numbers_OE

# Summarize events in CEx and IR + Alt3 and Alt5
EV_KD <- lapply(VTS_KD, function(x) Wrap_sepMIC(x)); Events=c("Alt3","Alt5","CEx","µEx","IR")
EV_OE <- lapply(VTS_OE, function(x) Wrap_sepMIC(x)); EV_KD; EV_OE
# EV_KD <- lapply(VTS_KD, function(x) Wrap_inclMIC(x)); Events=c("Alt3","Alt5","CEx","IR")
# EV_OE <- lapply(VTS_OE, function(x) Wrap_inclMIC(x)); 


## DIVIDE EVENTS IN UP / DOWN
UP_KD <- lapply(VTS_KD, function(x) {y <- subset(x, dPSI >0); return(y)})
UP_OE <- lapply(VTS_OE, function(x) {y <- subset(x, dPSI >0); return(y)})

DOWN_KD <- lapply(VTS_KD, function(x) {y <- subset(x, dPSI <0); return(y)})
DOWN_OE <- lapply(VTS_OE, function(x) {y <- subset(x, dPSI <0); return(y)})


## CALCULATE RELEVANT OVERLAPS/DIFFERENCES
# PROPORTION OF EVENTS THAT ARE OVERLAPPING IN THE 2 shRNAs
sum(is.element(VTS_KD$d12.d00_shC5$EVENT,VTS_KD$d12.d00_shC1$EVENT)) / length(VTS_KD$d12.d00_shC5$EVENT)
sum(is.element(VTS_KD$d12.d00_shU4$EVENT,VTS_KD$d12.d00_shU3$EVENT)) / length(VTS_KD$d12.d00_shU4$EVENT)
sum(is.element(VTS_KD$d12_shC1.shSCR$EVENT,VTS_KD$d12_shC5.shSCR$EVENT)) / length(VTS_KD$d12_shC1.shSCR$EVENT)
sum(is.element(VTS_KD$d12_shU3.shSCR$EVENT,VTS_KD$d12_shU4.shSCR$EVENT)) / length(VTS_KD$d12_shU3.shSCR$EVENT)

# PROPORTION OF EVENTS THAT CHANGE IN THE SAME DIRECTION IN THE 2 shRNAs CPSF3
sum(is.element(UP_KD$d12.d00_shC1$EVENT,UP_KD$d12.d00_shC5$EVENT)) / length(UP_KD$d12.d00_shC1$EVENT)
sum(is.element(DOWN_KD$d12.d00_shC1$EVENT,DOWN_KD$d12.d00_shC5$EVENT)) / length(DOWN_KD$d12.d00_shC1$EVENT)
(sum(is.element(UP_KD$d12.d00_shC1$EVENT,UP_KD$d12.d00_shC5$EVENT)) + sum(is.element(DOWN_KD$d12.d00_shC1$EVENT,DOWN_KD$d12.d00_shC5$EVENT))) / sum(is.element(VTS_KD$d12.d00_shC1$EVENT,VTS_KD$d12.d00_shC5$EVENT))
# PROPORTION OF EVENTS THAT CHANGE IN THE SAME DIRECTION IN THE 2 shRNAs UL1
sum(is.element(UP_KD$d12.d00_shU3$EVENT,UP_KD$d12.d00_shU4$EVENT)) / length(UP_KD$d12.d00_shU3$EVENT)
sum(is.element(DOWN_KD$d12.d00_shU3$EVENT,DOWN_KD$d12.d00_shU4$EVENT)) / length(DOWN_KD$d12.d00_shU3$EVENT)
(sum(is.element(UP_KD$d12.d00_shU3$EVENT,UP_KD$d12.d00_shU4$EVENT)) + sum(is.element(DOWN_KD$d12.d00_shU3$EVENT,DOWN_KD$d12.d00_shU4$EVENT))) / sum(is.element(VTS_KD$d12.d00_shU3$EVENT,VTS_KD$d12.d00_shU4$EVENT)) 

## RELEVANT EVENTS
EV_OVL_KD <- list(
  d12.d00_shCPSF3= VTS_KD$d12.d00_shC1$EVENT[is.element(VTS_KD$d12.d00_shC1$EVENT,VTS_KD$d12.d00_shC5$EVENT)],
  d12.d00_shUL1= VTS_KD$d12.d00_shU3$EVENT[is.element(VTS_KD$d12.d00_shU3$EVENT,VTS_KD$d12.d00_shU4$EVENT)]
)
# no # EV_OVL_KD$CPSF3_dep <-  EV_OVL_KD$d12.d00_shCPSF3[!is.element(EV_OVL_KD$d12.d00_shCPSF3,VTS_KD$d12.d00_shSCR$EVENT)]
# no # EV_OVL_KD$UL1_dep  <-  EV_OVL_KD$d12.d00_shUL1[!is.element(EV_OVL_KD$d12.d00_shUL1,VTS_KD$d12.d00_shSCR$EVENT)]
EV_OVL_KD$CPSF3_dep_d12 <-  VTS_KD$d12_shC1.shSCR$EVENT[is.element(VTS_KD$d12_shC1.shSCR$EVENT,VTS_KD$d12_shC5.shSCR$EVENT)]
EV_OVL_KD$UL1_dep_d12 <-  VTS_KD$d12_shU3.shSCR$EVENT[is.element(VTS_KD$d12_shU3.shSCR$EVENT,VTS_KD$d12_shU4.shSCR$EVENT)]

EV_OVL_OE <- list(
  d12.d00_TIA1= VTS_OE$d12.d00_TIA1$EVENT[is.element(VTS_OE$d12.d00_TIA1$EVENT,VTS_OE$d12.d00_TIA1$EVENT)]
)
# no # EV_OVL_OE$TIA1_dep = EV_OVL_OE$d12.d00_TIA1[!is.element(EV_OVL_OE$d12.d00_TIA1,VTS_OE$d12.d00_Empty$EVENT)]
EV_OVL_OE$TIA1_dep_d12= VTS_OE$d12_TIA1.Empty$EVENT

summary(EV_OVL_KD)
head(EV_OVL_KD$CPSF3_dep)
summary(EV_OVL_OE)
head(EV_OVL_OE$TIA1_dep)

## EXTRACT dPSI VALUES
dPSI_ctrls <- lapply(CTRLs, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")])
dPSI_KD <- lapply(VTS_KD, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")])
dPSI_OE <- lapply(VTS_OE, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")])
dPSI_KD_dPSI0 <- lapply(VTS_KD_dPSI0, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")])
dPSI_OE_dPSI0 <- lapply(VTS_OE_dPSI0, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")])



## PLOT dPSI VALUES KD
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/KD/Plots/")
shs <- c("d12.d00_shC1","d12.d00_shC5","d12.d00_shU3","d12.d00_shU4")
mdf_dPSI_KD <- melt(dPSI_KD[shs],level=1)
mdf_dPSI_KD_dPSI0 <- melt(dPSI_KD_dPSI0[shs],level=1)
colnames(mdf_dPSI_KD) = colnames(mdf_dPSI_KD_dPSI0) = c("EVENT","GENE","COMPLEX","variable","dPSI","Conditions")

# ## dPSI >= 10 upon any shRNA and dPSI >= 10 upon shSCR, but with a ddPSI of at least 10 
# mdf_dPSI_KD_SCRdPSI10 <- merge(x=mdf_dPSI_KD,y=dPSI_KD$d12.d00_shSCR,by="EVENT",suffixes = c("_test","_control"))
# mdf_dPSI_KD_SCRdPSI10$ddPSI <- mdf_dPSI_KD_SCRdPSI10$dPSI_test - mdf_dPSI_KD_SCRdPSI10$dPSI_control
# mdf_dPSI_KD_SCRdPSI10$abs_ddPSI <- abs(mdf_dPSI_KD_SCRdPSI10$ddPSI)
# mdf_dPSI_KD_SCRdPSI10[which(mdf_dPSI_KD_SCRdPSI10$abs_ddPSI <10),"abs_ddPSI"] <- NA
# head(mdf_dPSI_KD_SCRdPSI10)
# dim(mdf_dPSI_KD_SCRdPSI10)
# # Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_SCRdPSI10,name_x = "shSCR",name_y = "shRNA",suffix="",labels_ddPSI = 40)
# # Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_SCRdPSI10,name_x = "shSCR",name_y = "shRNA",suffix="nolabels",labels_ddPSI = 100)

# ## dPSI >= 10 upon any shRNA and no needed difference upon shSCR, but with a ddPSI of at least 10
# mdf_dPSI_KD_allSCR <- merge(x=mdf_dPSI_KD,y=dPSI_ctrls$d12.d00_shSCR_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
# mdf_dPSI_KD_allSCR$ddPSI <- mdf_dPSI_KD_allSCR$dPSI_test - mdf_dPSI_KD_allSCR$dPSI_control
# mdf_dPSI_KD_allSCR$abs_ddPSI <- abs(mdf_dPSI_KD_allSCR$ddPSI)
# mdf_dPSI_KD_allSCR[which(mdf_dPSI_KD_allSCR$abs_ddPSI <10),"abs_ddPSI"] <- NA
# head(mdf_dPSI_KD_allSCR)
# dim(mdf_dPSI_KD_allSCR)
# # Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_allSCR,name_x = "shSCR",name_y = "shRNA",suffix="allctrl",labels_ddPSI = 40)
# # Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_allSCR,name_x = "shSCR",name_y = "shRNA",suffix="allctrl_nolabels",labels_ddPSI = 100)

## dPSI >= 10 upon any shRNA and no needed difference upon shSCR, but with a ddPSI of at least 10
mdf_dPSI_KD_ovlshCU <- subset(mdf_dPSI_KD, EVENT %in% EV_OVL_KD$d12.d00_shCPSF3 | EVENT %in% EV_OVL_KD$d12.d00_shUL1 | EVENT %in% VTS_KD$d12.d00_shSCR)
mdf_dPSI_KD_allSCR_ovlshCU <- merge(x=mdf_dPSI_KD_ovlshCU,y=dPSI_ctrls$d12.d00_shSCR_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_KD_allSCR_ovlshCU$ddPSI <- mdf_dPSI_KD_allSCR_ovlshCU$dPSI_test - mdf_dPSI_KD_allSCR_ovlshCU$dPSI_control
mdf_dPSI_KD_allSCR_ovlshCU$abs_ddPSI <- abs(mdf_dPSI_KD_allSCR_ovlshCU$ddPSI)
mdf_dPSI_KD_allSCR_ovlshCU[which(mdf_dPSI_KD_allSCR_ovlshCU$abs_ddPSI <10),"abs_ddPSI"] <- NA
head(mdf_dPSI_KD_allSCR_ovlshCU)
dim(mdf_dPSI_KD_allSCR_ovlshCU)
# Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_allSCR_ovlshCU,name_x = "shSCR",name_y = "shRNA",suffix="ovl2shRNAs",labels_ddPSI = 40)
# Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_allSCR_ovlshCU,name_x = "shSCR",name_y = "shRNA",suffix="ovl2shRNAs_nolabels",labels_ddPSI = 100)

## No needed difference upon shSCR/shC/U, but with a ddPSI of at least 10
mdf_dPSI_KD_ovlshCU_dPSI0 <- subset(mdf_dPSI_KD_dPSI0, EVENT %in% EV_OVL_KD$d12.d00_shCPSF3 | EVENT %in% EV_OVL_KD$d12.d00_shUL1 | EVENT %in% VTS_KD$d12.d00_shSCR)
mdf_dPSI_KD_ovlshCU_dPSI0$Target <- str_replace(mdf_dPSI_KD_ovlshCU_dPSI0$Conditions,pattern = ".$","")
mdf_dPSI_KD_ovlshCU_dPSI0_av <- aggregate(dPSI ~ EVENT+GENE+COMPLEX+Target, mdf_dPSI_KD_ovlshCU_dPSI0, mean)
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- merge(x=mdf_dPSI_KD_ovlshCU_dPSI0_av, y=mdf_dPSI_KD_ovlshCU_dPSI0,by=c("EVENT","GENE","COMPLEX","Target"),suffixes = c("","_srep"),all=T)
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- merge(x=mdf_dPSI_KD_allSCR_ovlshCU_dPSI0,y=dPSI_ctrls$d12.d00_shSCR_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$ddPSI <- mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$dPSI_srep - mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$dPSI_control
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$abs_ddPSI <- abs(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$ddPSI)

### DEFINE CPSF3/UL1-dependent events accordingly
EV_OVL_KD$CPSF3_dep <- unique(subset(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0, abs_ddPSI >= 10 & grepl("d12.d00_shC",mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$Conditions))$EVENT)
EV_OVL_KD$UL1_dep <- unique(subset(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0, abs_ddPSI >= 10 & grepl("d12.d00_shU",mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$Conditions))$EVENT)
tmp <- subset(mdf_dPSI_KD_allSCR_ovlshCU,abs(dPSI_test) >=10 | abs(dPSI_control) >=10)
tmpc <- subset(tmp, grepl("d12.d00_shC",tmp$Conditions))
tmpu <- subset(tmp, grepl("d12.d00_shU",tmp$Conditions))
EV_OVL_KD$CPSF3_indep <- unique(subset(tmpc, abs(ddPSI) <= 2)$EVENT)
EV_OVL_KD$UL1_indep <- unique(subset(tmp, abs(ddPSI) <= 2)$EVENT)
lapply(EV_OVL_KD, function(x) length(x))

## ADD THESE CATEGORIES TO mdf
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$categC =mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$categU = rep(NA,nrow(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0))
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 %>%
  mutate(categC = replace(categC, EVENT %in% EV_OVL_KD$CPSF3_dep, "CPSF3_dep")) %>%
  mutate(categC = replace(categC, EVENT %in% EV_OVL_KD$CPSF3_indep, "CPSF3_indep"))  %>% 
  mutate(categU = replace(categU, EVENT %in% EV_OVL_KD$UL1_dep, "UL1_dep")) %>%
  mutate(categU = replace(categU, EVENT %in% EV_OVL_KD$UL1_indep, "UL1_indep")) 
head(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0)
dim(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0)

## PLOT SCATTERPLOT
setwd("~/Desktop/")
Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_allSCR_ovlshCU_dPSI0,name_x = "shSCR",name_y = "shRNA",suffix="ovl2shRNAa_nolabels",labels_ddPSI = 100)

## PLOT SCATTERPLOT onfly of CPSF3-UL1-dep events
ggplot(subset(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0, Target == "d12.d00_shC" & categC == "CPSF3_dep"),
       aes(x=dPSI_control,y=dPSI_test,color=abs_ddPSI)) +
  facet_wrap(~Conditions)+    # divide facets according to Conditions
  coord_equal(xlim = c(-100,100),ylim = c(-100,100)) +     # limits of x and y axes
  geom_vline(xintercept = c(-10,10),size=0.4,linetype = "dashed",color="gray60")+     # vertical lines (change it according to your dPSI_control threshold)
  geom_hline(yintercept = c(-10,10),size=0.4,linetype = "dashed",color="gray60")+     # horizontal lines (change it according to your dPSI_test threshold)
  geom_abline(slope = 1,intercept = c(-10,10), size=0.4,linetype = "solid",color="gray60")+     # diagonal lines (change it according to your ddPSI threshold)
  geom_point(data=subset(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0, Target == "d12.d00_shC" & categC == "CPSF3_indep"),size=1,color="grey80") +
  geom_point(size=1) +
  #geom_text_repel(data=labels,label=labels$EVENT, size=2) +    # labels of events
  scale_color_viridis(option='magma',direction = 1,begin=0.05) +     # color scale
  scale_y_continuous(name="shRNA") +
  scale_x_continuous(name="shSCR") +
  theme_classic() + 
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(),
         panel.grid.major.y=element_blank(),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
ggsave(filename = "Scatterplot_dPSI_CPSF3dep_indep.pdf",width=7,height=6,device = cairo_pdf)

ggplot(subset(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0, Target == "d12.d00_shU" & categU == "UL1_dep"),
       aes(x=dPSI_control,y=dPSI_test,color=abs_ddPSI)) +
  facet_wrap(~Conditions)+    # divide facets according to Conditions
  coord_equal(xlim = c(-100,100),ylim = c(-100,100)) +     # limits of x and y axes
  geom_vline(xintercept = c(-10,10),size=0.4,linetype = "dashed",color="gray60")+     # vertical lines (change it according to your dPSI_control threshold)
  geom_hline(yintercept = c(-10,10),size=0.4,linetype = "dashed",color="gray60")+     # horizontal lines (change it according to your dPSI_test threshold)
  geom_abline(slope = 1,intercept = c(-10,10), size=0.4,linetype = "solid",color="gray60")+     # diagonal lines (change it according to your ddPSI threshold)
  geom_point(data=subset(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0, Target == "d12.d00_shU" & categU == "UL1_indep"),size=1,color="grey80") +
  geom_point(size=1) +
  #geom_text_repel(data=labels,label=labels$EVENT, size=2) +    # labels of events
  scale_color_viridis(option='magma',direction = 1,begin=0.05) +     # color scale
  scale_y_continuous(name="shSNA") +
  scale_x_continuous(name="shSCR") +
  theme_classic() + 
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(),
         panel.grid.major.y=element_blank(),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
ggsave(filename = "Scatterplot_dPSI_UL1dep_indep.pdf",width=7,height=6,device = cairo_pdf)



## PLOT dPSI VALUES OE
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/OE/Plots/")
mdf_dPSI_OE <- melt(dPSI_OE["d12.d00_TIA1"],level=1)
mdf_dPSI_OE_dPSI0 <- melt(dPSI_OE_dPSI0["d12.d00_TIA1"],level=1)
colnames(mdf_dPSI_OE) = colnames(mdf_dPSI_OE_dPSI0) =c("EVENT","GENE","COMPLEX","variable","dPSI","Conditions")

# ## dPSI >= 10 upon TIA1 OE and dPSI >= 10 upon Empty, but with a ddPSI of at least 10 
# mdf_dPSI_OE_empdPSI10 <- merge(x=mdf_dPSI_OE,y=dPSI_OE$d12.d00_Empty,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
# mdf_dPSI_OE_empdPSI10$ddPSI <- mdf_dPSI_OE_empdPSI10$dPSI_test - mdf_dPSI_OE_empdPSI10$dPSI_control
# mdf_dPSI_OE_empdPSI10$abs_ddPSI <- abs(mdf_dPSI_OE_empdPSI10$ddPSI)
# mdf_dPSI_OE_empdPSI10[which(mdf_dPSI_OE_empdPSI10$abs_ddPSI <10),"abs_ddPSI"] <- NA
# head(mdf_dPSI_OE_empdPSI10)
# dim(mdf_dPSI_OE_empdPSI10)
# # Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_OE_empdPSI10,name_x = "Empty",name_y = "T7TIA1",suffix="",labels_ddPSI = 40)

## dPSI >= 10 upon TIA1 OE and no needed difference upon Empty, but with a ddPSI of at least 10
mdf_dPSI_OE_allEmp <- merge(x=mdf_dPSI_OE,y=dPSI_ctrls$d12.d00_Empty_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_OE_allEmp$ddPSI <- mdf_dPSI_OE_allEmp$dPSI_test - mdf_dPSI_OE_allEmp$dPSI_control
mdf_dPSI_OE_allEmp$abs_ddPSI <- abs(mdf_dPSI_OE_allEmp$ddPSI)
mdf_dPSI_OE_allEmp[which(mdf_dPSI_OE_allEmp$abs_ddPSI <10),"abs_ddPSI"] <- NA
head(mdf_dPSI_OE_allEmp)
dim(mdf_dPSI_OE_allEmp)
# Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_OE_allEmp,name_x = "Empty",name_y = "T7TIA1",suffix="allctrl",labels_ddPSI = 40)
# Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_OE_allEmp,name_x = "Empty",name_y = "T7TIA1",suffix="allctrl_nolabels",labels_ddPSI = 100)

## No needed difference upon Empty nor TIA1, but with a ddPSI of at least 10
mdf_dPSI_OE_allEmp_dPSI0 <- merge(x=mdf_dPSI_OE_dPSI0,y=dPSI_ctrls$d12.d00_Empty_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_OE_allEmp_dPSI0$ddPSI <- mdf_dPSI_OE_allEmp_dPSI0$dPSI_test - mdf_dPSI_OE_allEmp_dPSI0$dPSI_control
mdf_dPSI_OE_allEmp_dPSI0$abs_ddPSI <- abs(mdf_dPSI_OE_allEmp_dPSI0$ddPSI)
mdf_dPSI_OE_allEmp_dPSI0[which(mdf_dPSI_OE_allEmp_dPSI0$abs_ddPSI <10),"abs_ddPSI"] <- NA

head(mdf_dPSI_OE_allEmp_dPSI0)
dim(mdf_dPSI_OE_allEmp_dPSI0)

### DEFINE TIA1-dependent events accordingly
EV_OVL_OE$TIA1_dep <- unique(subset(mdf_dPSI_OE_allEmp_dPSI0, abs_ddPSI >= 10)$EVENT)
tmp <- subset(mdf_dPSI_OE_allEmp, abs(dPSI_test) >=10 | abs(dPSI_control) >=10)
EV_OVL_OE$TIA1_indep <- unique(subset(tmp, abs(ddPSI) <= 2)$EVENT)
lapply(EV_OVL_OE,function(x) length(x))
lapply(EV_OVL_OE,function(x) Wrap_sepMIC(subset(INCL$OE, EVENT %in% x)))

## ADD THESE CATEGORIES TO mdf
mdf_dPSI_OE_allEmp_dPSI0$categT  <-  rep(NA,nrow(mdf_dPSI_OE_allEmp_dPSI0))
mdf_dPSI_OE_allEmp_dPSI0 <- mdf_dPSI_OE_allEmp_dPSI0 %>%
  mutate(categT = replace(categT, EVENT %in% EV_OVL_OE$TIA1_dep, "TIA1_dep")) %>%
  mutate(categT = replace(categT, EVENT %in% EV_OVL_OE$TIA1_indep, "TIA1_indep"))
head(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0)
dim(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0)

## PLOT SCATTERPLOT
Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_OE_allEmp_dPSI0,name_x = "Empty",name_y = "T7TIA1",suffix="allctrl_nolabels_dPSI0",labels_ddPSI = 100)


## PLOT SCATTERPLOT only of TIA1-dep/indep events

ggplot(subset(mdf_dPSI_OE_allEmp_dPSI0, categT == "TIA1_dep"),
       aes(x=dPSI_control,y=dPSI_test,color=abs_ddPSI)) +
  #facet_wrap(~Target)+    # divide facets according to Conditions
  coord_equal(xlim = c(-100,100),ylim = c(-100,100)) +     # limits of x and y axes
  geom_vline(xintercept = c(-10,10),size=0.4,linetype = "dashed",color="gray60")+     # vertical lines (change it according to your dPSI_control threshold)
  geom_hline(yintercept = c(-10,10),size=0.4,linetype = "dashed",color="gray60")+     # horizontal lines (change it according to your dPSI_test threshold)
  geom_abline(slope = 1,intercept = c(-10,10), size=0.4,linetype = "solid",color="gray60")+     # diagonal lines (change it according to your ddPSI threshold)
  geom_point(data=subset(mdf_dPSI_OE_allEmp_dPSI0, categT == "TIA1_indep"),size=1,color="grey80") +
  geom_point(size=1) +
  #geom_text_repel(data=labels,label=labels$EVENT, size=2) +    # labels of events
  scale_color_viridis(option='magma',direction = 1,begin=0.05) +     # color scale
  scale_y_continuous(name="T7TIA1") +
  scale_x_continuous(name="Empty") +
  theme_classic() + 
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(),
         panel.grid.major.y=element_blank(),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
ggsave(filename = "Scatterplot_dPSI_TIA1dep_indep.pdf",width=7,height=6,device = cairo_pdf)




