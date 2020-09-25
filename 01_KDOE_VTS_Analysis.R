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
  d12.d00_shU4=read.table("KD/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t")
)
VTS_KD <- list(
  d12.d00_shSCR=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shSCR-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shC1=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n1-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shC5=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shCPSF3_n5-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shU3=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n3-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_shU4=read.table("KD/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-vs-KD_day12_shUL1_n4-with_dPSI.tab",header=T,sep="\t")
)
VTS_KD_dPSI0 <- lapply(VTS_KD_dPSI0, function(x) {rownames(x)=x$EVENT; x})  
VTS_KD <- lapply(VTS_KD, function(x) {rownames(x)=x$EVENT; x})  
lapply(VTS_KD_dPSI0, function(x) dim(x))  
lapply(VTS_KD, function(x) dim(x))  

## IMPORT DIFFAS DATA FOR OE
VTS_OE_dPSI0 <- list(
  d12.d00_Empty=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_TIA1=read.table("OE/DiffAS_dPSI0/DiffAS-Mm218-dPSI0-range0-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t")
)
VTS_OE <- list(
  d12.d00_Empty=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_Empty-with_dPSI.tab",header=T,sep="\t"),
  d12.d00_TIA1=read.table("OE/DiffAS_dPSI10/DiffAS-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-vs-OE_day12_T7TIA1-with_dPSI.tab",header=T,sep="\t")
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
EV_KD <- lapply(VTS_KD, function(x) Wrap_VTS_Events(x, "out_CEx")); #Events=c("Alt3","Alt5","IR","CEx","µEx")
EV_OE <- lapply(VTS_OE, function(x) Wrap_VTS_Events(x, "out_CEx")); EV_KD; EV_OE


## DIVIDE EVENTS IN UP / DOWN
UP_KD <- lapply(VTS_KD, function(x) {y <- subset(x, dPSI >0); return(y)})
UP_OE <- lapply(VTS_OE, function(x) {y <- subset(x, dPSI >0); return(y)})

DOWN_KD <- lapply(VTS_KD, function(x) {y <- subset(x, dPSI <0); return(y)})
DOWN_OE <- lapply(VTS_OE, function(x) {y <- subset(x, dPSI <0); return(y)})


## CALCULATE RELEVANT OVERLAPS/DIFFERENCES
# PROPORTION OF EVENTS THAT ARE OVERLAPPING IN THE 2 shRNAs
sum(is.element(VTS_KD$d12.d00_shC1$EVENT,VTS_KD$d12.d00_shC5$EVENT)) / length(VTS_KD$d12.d00_shC1$EVENT)
sum(is.element(VTS_KD$d12.d00_shU3$EVENT,VTS_KD$d12.d00_shU4$EVENT)) / length(VTS_KD$d12.d00_shU3$EVENT)

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

EV_OVL_OE <- list(
  d12.d00_TIA1= VTS_OE$d12.d00_TIA1$EVENT[is.element(VTS_OE$d12.d00_TIA1$EVENT,VTS_OE$d12.d00_TIA1$EVENT)]
)

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

## No needed dPSI upon shSCR/shC/U, but with a ddPSI of at least 10
mdf_dPSI_KD_ovlshCU_dPSI0 <- subset(mdf_dPSI_KD_dPSI0, EVENT %in% EV_OVL_KD$d12.d00_shCPSF3 | EVENT %in% EV_OVL_KD$d12.d00_shUL1 | EVENT %in% VTS_KD$d12.d00_shSCR)
mdf_dPSI_KD_ovlshCU_dPSI0$Target <- str_replace(mdf_dPSI_KD_ovlshCU_dPSI0$Conditions,pattern = ".$","")
mdf_dPSI_KD_ovlshCU_dPSI0_av <- aggregate(dPSI ~ EVENT+GENE+COMPLEX+Target, mdf_dPSI_KD_ovlshCU_dPSI0, mean)
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- merge(x=mdf_dPSI_KD_ovlshCU_dPSI0_av, y=mdf_dPSI_KD_ovlshCU_dPSI0,by=c("EVENT","GENE","COMPLEX","Target"),suffixes = c("","_srep"),all=T)
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0 <- merge(x=mdf_dPSI_KD_allSCR_ovlshCU_dPSI0,y=dPSI_ctrls$d12.d00_shSCR_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$ddPSI <- mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$dPSI_srep - mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$dPSI_control
mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$abs_ddPSI <- abs(mdf_dPSI_KD_allSCR_ovlshCU_dPSI0$ddPSI)

## dPSI >= 10 upon any shRNA and no needed difference upon shSCR, but with a ddPSI of at least 10
mdf_dPSI_KD_ovlshCU <- subset(mdf_dPSI_KD, EVENT %in% EV_OVL_KD$d12.d00_shCPSF3 | EVENT %in% EV_OVL_KD$d12.d00_shUL1 | EVENT %in% VTS_KD$d12.d00_shSCR)
mdf_dPSI_KD_allSCR_ovlshCU <- merge(x=mdf_dPSI_KD_ovlshCU,y=dPSI_ctrls$d12.d00_shSCR_all,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_KD_allSCR_ovlshCU$ddPSI <- mdf_dPSI_KD_allSCR_ovlshCU$dPSI_test - mdf_dPSI_KD_allSCR_ovlshCU$dPSI_control
mdf_dPSI_KD_allSCR_ovlshCU$abs_ddPSI <- abs(mdf_dPSI_KD_allSCR_ovlshCU$ddPSI)
mdf_dPSI_KD_allSCR_ovlshCU[which(mdf_dPSI_KD_allSCR_ovlshCU$abs_ddPSI <10),"abs_ddPSI"] <- NA
head(mdf_dPSI_KD_allSCR_ovlshCU)
dim(mdf_dPSI_KD_allSCR_ovlshCU)

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
# Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_KD_allSCR_ovlshCU_dPSI0,name_x = "shSCR",name_y = "shRNA",suffix="ovl2shRNAa_nolabels",labels_ddPSI = 100)

## PLOT SCATTERPLOT only of CPSF3-UL1-dep events
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
# ggsave(filename = "Scatterplot_dPSI_CPSF3dep_indep.pdf",width=7,height=6,device = cairo_pdf)

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
# ggsave(filename = "Scatterplot_dPSI_UL1dep_indep.pdf",width=7,height=6,device = cairo_pdf)



## PLOT dPSI VALUES OE
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/OE/Plots/")
mdf_dPSI_OE <- melt(dPSI_OE["d12.d00_TIA1"],level=1)
mdf_dPSI_OE_dPSI0 <- melt(dPSI_OE_dPSI0["d12.d00_TIA1"],level=1)
colnames(mdf_dPSI_OE) = colnames(mdf_dPSI_OE_dPSI0) =c("EVENT","GENE","COMPLEX","variable","dPSI","Conditions")

## No needed dPSI upon Empty nor TIA1, but with a ddPSI of at least 10
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
# Plot_ScatterPlots_dPSI(mdf_dPSI = mdf_dPSI_OE_allEmp_dPSI0,name_x = "Empty",name_y = "T7TIA1",suffix="allctrl_nolabels_dPSI0",labels_ddPSI = 100)


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
# ggsave(filename = "Scatterplot_dPSI_TIA1dep_indep.pdf",width=7,height=6,device = cairo_pdf)




## GENERATE TABLE WITH SETS OF EVENTS FOR KD
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/KD/")
Final_KD <- lapply(EV_OVL_KD, function(x) subset(INCL$KD,EVENT %in% x) )
# Add SCR control events
Final_KD$d12.d00_shSCR <- VTS_KD$d12.d00_shSCR
lapply(Final_KD,function(x) dim(x))
## EXPORT TABLES
# for (n in names(Final_KD)) {
#   write.table(Final_KD[[n]], file=paste("INCL_",n,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = NA)
# }

## GENERATE TABLE WITH SETS OF EVENTS FOR OE
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/OE/")
Final_OE <- lapply(EV_OVL_OE, function(x) subset(INCL$OE,EVENT %in% x) )
# Add Empty control events
Final_OE$d12.d00_Empty <- VTS_OE$d12.d00_Empty
lapply(Final_OE,function(x) dim(x))
## EXPORT TABLES
# for (n in names(Final_OE)) {
# write.table(Final_OE[[n]], file=paste("INCL_",n,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = NA)
# }


## NUMBER OF EVENTS
lapply(Final_KD, function(x) Wrap_VTS_Events(x, "in_CEx"))
lapply(Final_OE, function(x) Wrap_VTS_Events(x, "in_CEx"))



## OVERLAP CPSF3/UL1/TIA1-dep events
overlapCUT <- calculate.overlap(list(EV_OVL_KD$CPSF3_dep,EV_OVL_KD$UL1_dep,EV_OVL_OE$TIA1_dep))
# pdf(file="VennDiagram_CPSF3_UL1_TIA1dep.pdf",width = 8, height=8)
draw.triple.venn( area1 = length(overlapCUT$a1),
                  area2 = length(overlapCUT$a2),
                  area3 = length(overlapCUT$a3),
                  n12 = length(overlapCUT$a2),
                  n23 = length(overlapCUT$a2),
                  n13 = length(overlapCUT$a2),
                  n123 = length(overlapCUT$a5),
                  cross.area = length(overlapCU$a3),
                  category = c("CPSF3-dependent","UL1-dependent"),
                  lwd=0,
                  fill=c("darkslategray","maroon"),
                  cat.fontface = "bold", cex = 3,
                  cat.fontfamily = "Helvetica",
                  cat.pos =0,
                  #print.mode=c('raw','percent'),
                  fontfamily = "Helvetica"
)
# dev.off()

## OVERLAP CPSF3-UL1 in same direction?
sh_C <- c("d12.d00_shC1","d12.d00_shC5")
sh_U <- c("d12.d00_shU3","d12.d00_shU4")
mdf_dPSI_KD_C <- melt(dPSI_KD[sh_C],level=1,value.name="dPSI")
  mdf_dPSI_KD_C_av <- aggregate(dPSI ~ EVENT+GENE+COMPLEX, mdf_dPSI_KD_C, mean)
mdf_dPSI_KD_U <- melt(dPSI_KD[sh_U],level=1,value.name="dPSI")
  mdf_dPSI_KD_U_av <- aggregate(dPSI ~ EVENT+GENE+COMPLEX, mdf_dPSI_KD_U, mean)

mdf_dPSI_KD_CU <- merge(x=mdf_dPSI_KD_C_av,y=mdf_dPSI_KD_U_av,by=c("EVENT","GENE","COMPLEX"),suffixes = c("_test","_control"))
mdf_dPSI_KD_CU$ddPSI <- mdf_dPSI_KD_CU$dPSI_test - mdf_dPSI_KD_CU$dPSI_control
mdf_dPSI_KD_CU$Conditions <- rep("",nrow(mdf_dPSI_KD_CU))
mdf_dPSI_KD_CU$abs_ddPSI <- abs(mdf_dPSI_KD_CU$ddPSI)
head(mdf_dPSI_KD_CU)
dim(mdf_dPSI_KD_CU)

# How many go in the same direction (99.6%)
sum(sign(mdf_dPSI_KD_CU$dPSI_test)*sign(mdf_dPSI_KD_CU$dPSI_control)>0)
# How many have opposite sign (0.03%)
sum(sign(mdf_dPSI_KD_CU$dPSI_test)*sign(mdf_dPSI_KD_CU$dPSI_control)<0)



## PLOT VIOLIN PLOTs WITH PSI VALUES OF ALL KD SETS
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/KD/Plots/")
sep_Final_KD <- SplitList_inclMIC(Final_KD)
for(n in names(sep_Final_KD)) {
  tmp <- sep_Final_KD[[n]]
  df_KD <- lapply(tmp, function(x) x[,!grepl(".Q",colnames(x))])
  df_KD <- lapply(df_KD, function(x) x[,!grepl("^dPSI",colnames(x))])
  df_KD <- lapply(df_KD, function(x) x[,!grepl("^OE_",colnames(x))])
  df_KD$d12.d00_shCPSF3 <- NULL; df_KD$d12.d00_shUL1 <- NULL
  mdf_PSIs_KD <- melt(df_KD, id.vars =c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX"), 
                      variable.name = "Sample", value.name = "PSI",
                      level=1) 
  # }
  mdf_PSIs_KD$PSI <- as.numeric(mdf_PSIs_KD$PSI)
  mdf_PSIs_KD$Condition <- str_replace(pattern = "_[12]$",replacement = "",string = mdf_PSIs_KD$Sample)
  mdf_PSIs_KD <- aggregate (PSI ~ GENE+EVENT+COORD+LENGTH+FullCO+COMPLEX+L1+Condition, mdf_PSIs_KD, mean)
  mdf_PSIs_KD$Condition <- factor(mdf_PSIs_KD$Condition, levels = c("KD_day00_Ctrl","KD_day12_shSCR","KD_day12_shCPSF3_n1","KD_day12_shCPSF3_n5","KD_day12_shUL1_n3","KD_day12_shUL1_n4"))
  mdf_PSIs_KD$categ = rep(NA,nrow(mdf_PSIs_KD))  
  mdf_PSIs_KD <- mdf_PSIs_KD %>%
    mutate(categ = replace(categ, PSI > 75 | PSI < 25, "non-AS")) %>%
    mutate(categ = replace(categ, PSI >= 25 & PSI <= 75, "AS"))
  head(mdf_PSIs_KD)
  
  # for(u in unique(mdf_PSIs_KD$L1)) {
  #   x <- subset(mdf_PSIs_KD, L1 == u)
  # ggplot(x,
  #        aes(x=Condition, y=PSI, fill=Condition)) +
  #   geom_violin(trim = T,color="transparent",show.legend = F) +
  #   geom_boxplot(width=0.1,color="grey20",show.legend = F,outlier.size=0.3) +
  #   scale_x_discrete() +
  #   scale_fill_viridis(option = "D",discrete = T,direction = -1,alpha = 0.55) +
  #     ggtitle(paste0(n," _ ",u,", n=",length(unique(x$EVENT)))) +
  #     theme_classic() +
  #   theme (text = element_text(color="grey20",size=10),
  #          axis.title = element_text(face="bold"),
  #          axis.text.x = element_text(angle = 30,hjust=1,vjust=1,size =7),
  #          legend.position = "bottom", legend.title = element_blank(),
  #          panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
  #          strip.background = element_blank(),strip.text = element_text(face="bold"))
  # ggsave(filename = paste0("Violin_",n,"_",u,"_av.pdf"),width=7,height=6,device = cairo_pdf)
  # }
}

## PLOT VIOLIN PLOTs WITH PSI VALUES OF OE SETS
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/OE/Plots/")
sep_Final_OE <- SplitList_inclMIC(Final_OE)
for (n in names(sep_Final_OE)) {
  tmp <- sep_Final_OE[[n]]
  df_OE <- lapply(tmp, function(x) x[,!grepl(".Q",colnames(x))])
  df_OE <- lapply(df_OE, function(x) x[,!grepl("^dPSI",colnames(x))])
  df_OE <- lapply(df_OE, function(x) x[,!grepl("^KD_",colnames(x))])
  df_OE$d12.d00_TIA1 <- NULL }
mdf_PSIs_OE <- melt(df_OE, id.vars =c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX"), 
                    variable.name = "Sample", value.name = "PSI",
                    level=1)
mdf_PSIs_OE$PSI <- as.numeric(mdf_PSIs_OE$PSI)
mdf_PSIs_OE$Condition <- str_replace(pattern = "_[12]$",replacement = "",string = mdf_PSIs_OE$Sample)
mdf_PSIs_OE <- aggregate (PSI ~ GENE+EVENT+COORD+LENGTH+FullCO+COMPLEX+L1+Condition, mdf_PSIs_OE, mean)
mdf_PSIs_OE$Condition <- factor(mdf_PSIs_OE$Condition, levels = c("OE_day00_Ctrl","OE_day12_Empty","OE_day12_T7TIA1"))
mdf_PSIs_OE$categ = rep(NA,nrow(mdf_PSIs_OE))  
mdf_PSIs_OE <- mdf_PSIs_OE %>%
  mutate(categ = replace(categ, PSI > 75 | PSI < 25, "non-AS")) %>%
  mutate(categ = replace(categ, PSI >= 25 & PSI <= 75, "AS"))
head(mdf_PSIs_OE)

# for(u in unique(mdf_PSIs_OE$L1)) {
#   x <- subset(mdf_PSIs_OE, L1 == u)
#   ggplot(x, 
#          aes(x=Condition, y=PSI, fill=Condition)) +
#     geom_violin(trim = T,color="transparent",show.legend = F) +
#     geom_boxplot(width=0.1,color="grey20",show.legend = F,outlier.size=0.3) +
#     scale_x_discrete() +
#     scale_fill_viridis(option = "D",discrete = T,direction = -1,alpha = 0.55) +
#     ggtitle(paste0(n," _ ",u,", n=",length(unique(x$EVENT)))) + 
#     theme_classic() +
#     theme (text = element_text(color="grey20",size=10),
#            axis.title = element_text(face="bold"), 
#            axis.text.x = element_text(angle = 30,hjust=1,vjust=1,size =7),
#            legend.position = "bottom", legend.title = element_blank(),
#            panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
#            strip.background = element_blank(),strip.text = element_text(face="bold"))
#   ggsave(filename = paste0("Violin_",n,"_",u,"_av.pdf"),width=7,height=6,device = cairo_pdf)
# }
#}



## OVERLAP WITH B CELL REPROGRAMMING AS EVENTS
## percentage of all events overlapping in general with B/M clusters
100*(sum(is.element(EV_OVL_KD$CPSF3_dep,rownames(B_PSIs_VTS_av$dPSI10))) )/length(EV_OVL_KD$CPSF3_dep)
100*(sum(is.element(EV_OVL_KD$UL1_dep,rownames(B_PSIs_VTS_av$dPSI10))) )/length(EV_OVL_KD$UL1_dep)
100*(sum(is.element(EV_OVL_OE$TIA1_dep,rownames(B_PSIs_VTS_av$dPSI10))) )/length(EV_OVL_OE$TIA1_dep)


ovl_allevents <- data.frame(row.names=c(names(Final_KD),names(Final_OE)))
for (e in names(Final_KD)) {
  ovl_allevents[e,"Condition"] = e
  ovl_allevents[e,"B_yes"] = 100*(sum(is.element(rownames(Final_KD[[e]]),rownames(B_PSIs_VTS_av$dPSI10))) )/nrow(Final_KD[[e]])
  ovl_allevents[e,"B_no"] = 100 - 100*(sum(is.element(rownames(Final_KD[[e]]),rownames(B_PSIs_VTS_av$dPSI10))) )/nrow(Final_KD[[e]])
  ovl_allevents[e,"M_yes"] = 100*(sum(is.element(rownames(Final_KD[[e]]),rownames(M_PSIs_VTS_av$dPSI10))) )/nrow(Final_KD[[e]])
  ovl_allevents[e,"M_no"] = 100 - 100*(sum(is.element(rownames(Final_KD[[e]]),rownames(M_PSIs_VTS_av$dPSI10))) )/nrow(Final_KD[[e]])
}
for (e in names(Final_OE)) {
  ovl_allevents[e,"Condition"] = e
  ovl_allevents[e,"B_yes"] = 100*(sum(is.element(rownames(Final_OE[[e]]),rownames(B_PSIs_VTS_av$dPSI10))) )/nrow(Final_OE[[e]])
  ovl_allevents[e,"B_no"] = 100 - 100*(sum(is.element(rownames(Final_OE[[e]]),rownames(B_PSIs_VTS_av$dPSI10))) )/nrow(Final_OE[[e]])
  ovl_allevents[e,"M_yes"] = 100*(sum(is.element(rownames(Final_OE[[e]]),rownames(M_PSIs_VTS_av$dPSI10))) )/nrow(Final_OE[[e]])
  ovl_allevents[e,"M_no"] = 100 - 100*(sum(is.element(rownames(Final_OE[[e]]),rownames(M_PSIs_VTS_av$dPSI10))) )/nrow(Final_OE[[e]])
}
ovl_allevents
mdf_ovl_allev <- melt(ovl_allevents, variable.name = "Overlap",value.name = "Percentage")
head(mdf_ovl_allev)

## PLOT BARPLOT WITH PERCENTAGES OF OVERLAP WITH B CELL REPROGRAMMING AS EVENTS
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/KD/Overlap_BMreprog/")
toplot <- subset(mdf_ovl_allev, 
                 grepl("^B_",mdf_ovl_allev$Overlap) & Condition %in% c("TIA1_dep","CPSF3_dep","UL1_dep","TIA1_indep","CPSF3_indep","UL1_indep","d12.d00_Empty","d12.d00_shSCR"))
toplot$Condition <- factor(toplot$Condition, levels= c("TIA1_dep","CPSF3_dep","UL1_dep","TIA1_indep","CPSF3_indep","UL1_indep","d12.d00_Empty","d12.d00_shSCR"))
p <- ggplot(toplot, aes (x=Condition,y= Percentage, fill = Overlap)) +
  geom_bar(width = 0.8,
           position=position_stack(reverse=T), 
           stat="identity") 
plot_tab = ggplot_build(p)$data[[1]]
plot_tab = plot_tab[order(plot_tab$group),]
plot_tab = plot_tab[order(plot_tab$fill,decreasing = T),]
toplot$coord <- apply(plot_tab, 1, function(x) mean(as.numeric(x[c("ymax","ymin")])))

p + 
  scale_fill_manual(values = c("#2e4e4e","grey80")) +    ##c("#2e4e4e","#af305f")
  scale_y_continuous(limits=c(0,100))+
  geom_text(data=toplot, aes(x= Condition, y=coord,
                             label=paste(round(Percentage,2),"%",sep=" ") ,
                             hjust=0.5,vjust=0.5,size=7),
            position = position_dodge(width=0),size=3) + 
  theme_minimal() +
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=45,hjust=1,vjust=1),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_line(linetype = "dashed",color="gray80"))
# ggsave(filename = paste0("Overlap_B_depindepAll.pdf"),width=7,height=6,device = cairo_pdf)
