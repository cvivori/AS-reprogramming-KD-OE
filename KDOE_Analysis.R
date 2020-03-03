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
VTS_KD <- list(
  dPSI00 = list(
      d12.d00_shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shSCR-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shC1=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shC5=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shU3=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shU4=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t"),
      d12_shC1.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
      d12_shC5.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
      d12_shU3.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
      d12_shU4.shSCR=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t")  ),
  dPSI10 = list(
      d12.d00_shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shSCR-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shC1=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shC5=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shU3=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_shU4=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t"),
      d12_shC1.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
      d12_shC5.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
      d12_shU3.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
      d12_shU4.shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day12_shSCR-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t") )
)
VTS_KD$dPSI00 <- lapply(VTS_KD$dPSI00, function(x) {rownames(x)=x$EVENT; x})  
VTS_KD$dPSI10 <- lapply(VTS_KD$dPSI10, function(x) {rownames(x)=x$EVENT; x})  
lapply(VTS_KD$dPSI00, function(x) dim(x))  
lapply(VTS_KD$dPSI10, function(x) dim(x))  


## IMPORT DIFFAS DATA FOR OE
VTS_OE <- list(
  dPSI00 = list(
      d12.d00_Empty=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_TIA1=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
      d12_TIA1.Empty=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day12_Empty-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t") ),
  dPSI10 = list(
      d12.d00_Empty=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_Empty-with_dPSI.tab",header=T,sep="\t"),
      d12.d00_TIA1=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
      d12_TIA1.Empty=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day12_Empty-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t") )
)
VTS_OE$dPSI00 <- lapply(VTS_OE$dPSI00, function(x) {rownames(x)=x$EVENT; x})
VTS_OE$dPSI10 <- lapply(VTS_OE$dPSI10, function(x) {rownames(x)=x$EVENT; x})
lapply(VTS_OE$dPSI00, function(x) dim(x))  
lapply(VTS_OE$dPSI10, function(x) dim(x))  



## TYPES OF EVENTS
# Find duplicate lines with FullCO  (Alt3 and Alt5)
Pairs=lapply(VTS_KD$dPSI10, function(x) FindPairs(table=x, column = "FullCO"))  
# Check if all duplicates are ALT3 or ALT3
lapply(Pairs, function(x) CheckAlt(x))

# Find duplicate lines with FullCO  (Alt3 and Alt5)
Pairs=lapply(VTS_OE$dPSI10, function(x) FindPairs(table=x, column = "FullCO"))  
# Check if all duplicates are ALT3 or ALT3
lapply(Pairs, function(x) CheckAlt(x))

# Calculate how many events in each table
Numbers_KD=lapply(VTS_KD$dPSI10, function(x) summary(x[,"COMPLEX"]))
Numbers_OE=lapply(VTS_OE$dPSI10, function(x) summary(x[,"COMPLEX"]))
Numbers_KD ; Numbers_OE

# Summarize events in CEx and IR + Alt3 and Alt5
EV_KD <- lapply(VTS_KD$dPSI10, function(x) Wrap_VTS_Events(x))
EV_OE <- lapply(VTS_OE$dPSI10, function(x) Wrap_VTS_Events(x)); EV_KD; EV_OE

## DIVIDE EVENTS IN UP / DOWN
UP_KD <- lapply(VTS_KD, function(y) lapply(y, function(x) {z <- subset(x, dPSI >0); return(z)}))
UP_OE <- lapply(VTS_OE, function(y) lapply(y, function(x) {z <- subset(x, dPSI >0); return(z)}))

DO_KD <- lapply(VTS_KD, function(y) lapply(y, function(x) {z <- subset(x, dPSI <0); return(z)}))
DO_OE <- lapply(VTS_OE, function(y) lapply(y, function(x) {z <- subset(x, dPSI <0); return(z)}))
                                                


## CALCULATE RELEVANT OVERLAPS/DIFFERENCES
# PROPORTION OF EVENTS THAT ARE OVERLAPPING IN THE 2 shRNAs
sum(is.element(VTS_KD$dPSI10$d12.d00_shC1$EVENT,VTS_KD$dPSI10$d12.d00_shC5$EVENT)) / length(VTS_KD$dPSI10$d12.d00_shC1$EVENT)
sum(is.element(VTS_KD$dPSI10$d12.d00_shU3$EVENT,VTS_KD$dPSI10$d12.d00_shU4$EVENT)) / length(VTS_KD$dPSI10$d12.d00_shU3$EVENT)
sum(is.element(VTS_KD$dPSI10$d12_shC1.shSCR$EVENT,VTS_KD$dPSI10$d12_shC5.shSCR$EVENT)) / length(VTS_KD$dPSI10$d12_shC1.shSCR$EVENT)
sum(is.element(VTS_KD$dPSI10$d12_shU3.shSCR$EVENT,VTS_KD$dPSI10$d12_shU4.shSCR$EVENT)) / length(VTS_KD$dPSI10$d12_shU3.shSCR$EVENT)

# PROPORTION OF EVENTS THAT CHANGE IN THE SAME DIRECTION IN THE 2 shRNAs CPSF3
sum(is.element(UP_KD$dPSI10$d12.d00_shC1$EVENT, UP_KD$dPSI10$d12.d00_shC5$EVENT)) / length(UP_KD$dPSI10$d12.d00_shC1$EVENT)
sum(is.element(DO_KD$dPSI10$d12.d00_shC1$EVENT, DO_KD$dPSI10$d12.d00_shC5$EVENT)) / length(DO_KD$dPSI10$d12.d00_shC1$EVENT)
# PROPORTION OF EVENTS THAT CHANGE IN THE SAME DIRECTION IN THE 2 shRNAs UL1
sum(is.element(UP_KD$dPSI10$d12.d00_shU3$EVENT, UP_KD$dPSI10$d12.d00_shU4$EVENT)) / length(UP_KD$dPSI10$d12.d00_shU3$EVENT)
sum(is.element(DO_KD$dPSI10$d12.d00_shU3$EVENT, DO_KD$dPSI10$d12.d00_shU4$EVENT)) / length(DO_KD$dPSI10$d12.d00_shU3$EVENT)



#############################################################################################################################

## SELECT EVENTS CHANGING day12-day0 IN BOTH shRNAs
EV_OVL_KD <- list()
EV_OVL_KD$d12.d00_shCPSF3= VTS_KD$dPSI10$d12.d00_shC1$EVENT[is.element(VTS_KD$dPSI10$d12.d00_shC1$EVENT, VTS_KD$dPSI10$d12.d00_shC5$EVENT)]
EV_OVL_KD$d12.d00_shUL1= VTS_KD$dPSI10$d12.d00_shU3$EVENT[is.element(VTS_KD$dPSI10$d12.d00_shU3$EVENT, VTS_KD$dPSI10$d12.d00_shU4$EVENT)]
## SELECT EVENTS CHANGING day12-day12shSCR IN BOTH shRNAs
EV_OVL_KD$CPSF3_dep_d12 <-  VTS_KD$dPSI10$d12_shC1.shSCR$EVENT[is.element(VTS_KD$dPSI10$d12_shC1.shSCR$EVENT, VTS_KD$dPSI10$d12_shC5.shSCR$EVENT)]
EV_OVL_KD$UL1_dep_d12 <-  VTS_KD$dPSI10$d12_shU3.shSCR$EVENT[is.element(VTS_KD$dPSI10$d12_shU3.shSCR$EVENT, VTS_KD$dPSI10$d12_shU4.shSCR$EVENT)]

## SELECT EVENTS CHANGING day12-day0 AND day12-day12Empty
EV_OVL_OE <- list()
EV_OVL_OE$d12.d00_TIA1= VTS_OE$dPSI10$d12.d00_TIA1$EVENT
EV_OVL_OE$TIA1_dep_d12= VTS_OE$dPSI10$d12_TIA1.Empty$EVENT


## EXTRACT dPSI VALUES
dPSI_ctrls <- lapply(CTRLs, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")])
dPSI_KD <- lapply(VTS_KD, function(y) lapply(y, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")]))
dPSI_OE <- lapply(VTS_OE, function(y) lapply(y, function(x) x[,c("EVENT","GENE","COMPLEX","dPSI")]))




## DEFINE CPSF3/UL1-dependent AND -independent EVENTS
shs <- c("d12.d00_shC1","d12.d00_shC5","d12.d00_shU3","d12.d00_shU4")
mdf_dPSI_KD <- lapply(dPSI_KD, function(x) melt(x[shs],level=1))
mdf_dPSI_KD <- lapply(mdf_dPSI_KD, function(x) {colnames(x) = c("EVENT","GENE","COMPLEX","variable","dPSI","Conditions"); x})
head(mdf_dPSI_KD$dPSI10)
lapply(mdf_dPSI_KD,dim)

## No needed difference upon shSCR/shC/U, but with a ddPSI of at least 10
mdf_dPSI_KD$dPSI00 %>% filter(EVENT %in% EV_OVL_KD$d12.d00_shCPSF3 | EVENT %in% EV_OVL_KD$d12.d00_shUL1)
dim( subset(mdf_dPSI_KD$dPSI00, EVENT %in% EV_OVL_KD$d12.d00_shCPSF3 | EVENT %in% EV_OVL_KD$d12.d00_shUL1))
mdf_dPSI_KD_ovlshCU_dPSI0$Target <- str_replace(mdf_dPSI_KD_ovlshCU_dPSI0$Conditions,pattern = ".$","")
mdf_dPSI_KD_ovlshCU_dPSI0_av <- aggregate(dPSI ~ EVENT+GENE+COMPLEX+Target, mdf_dPSI_KD_ovlshCU_dPSI0, mean)
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- merge(x=mdf_dPSI_KD_ovlshCU_dPSI0_av, y=mdf_dPSI_KD_ovlshCU_dPSI0,by=c("EVENT","GENE","COMPLEX","Target"),suffixes = c("","_srep"),all=T)
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- merge(x=mdf_dPSI_KD_allSCR_ovlshCU_dPSI0,y=dPSI_ctrls$d12.d00_shSCR_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$ddPSI <- mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$dPSI_srep - mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$dPSI_control
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$abs_ddPSI <- abs(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$ddPSI)


setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/KD/Plots/")




