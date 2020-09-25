require(tidyverse)
require(ggrepel)

## IMPORT
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/GEx_EdgeR/")
KDs_edgeR <- read.table("KD/FilteredGeneCpm_KD.txt", header=T, sep="\t")
OEs_edgeR <- read.table("OE/FilteredGeneCpm_OE.txt", header=T, sep="\t")
head(KDs_edgeR)


## CHECK KD EFFICIENCY
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/GEx_EdgeR/KD/")
KDeffic <- subset(KDs_edgeR, GeneName %in% c("Cpsf3","Hnrnpul1"))       
mKDeffic <- melt(KDeffic, variable.name = "Sample", value.name = "CPMs")
mKDeffic$Condition <- gsub("_1","",mKDeffic$Sample)
mKDeffic$Condition <- gsub("_2","",mKDeffic$Condition)
mKDeffic$Condition <- factor(mKDeffic$Condition, levels =c("KD_day00_Ctrl","KD_day12_shSCR","KD_day12_shCPSF3n1","KD_day12_shCPSF3n5","KD_day12_shUL1n3","KD_day12_shUL1n4"))
mKDeffic$Replicate <- str_sub(mKDeffic$Sample,-1)
mKDeffic$shRNA <- str_replace(pattern = "n[1345]", replacement = "", string = mKDeffic$Condition)
head(mKDeffic)

mKDeffic_av <- aggregate(CPMs ~ ENSGID+GeneName+shRNA, mKDeffic, mean)
subset(mKDeffic_av, GeneName == "Cpsf3" & shRNA == "KD_day12_shCPSF3")$CPMs / subset(mKDeffic_av, GeneName == "Cpsf3" & shRNA == "KD_day12_shSCR")$CPMs 
subset(mKDeffic_av, GeneName == "Cpsf3" & shRNA == "KD_day12_shUL1")$CPMs / subset(mKDeffic_av, GeneName == "Cpsf3" & shRNA == "KD_day12_shSCR")$CPMs 
subset(mKDeffic_av, GeneName == "Hnrnpul1" & shRNA == "KD_day12_shCPSF3")$CPMs / subset(mKDeffic_av, GeneName == "Hnrnpul1" & shRNA == "KD_day12_shSCR")$CPMs 
subset(mKDeffic_av, GeneName == "Hnrnpul1" & shRNA == "KD_day12_shUL1")$CPMs / subset(mKDeffic_av, GeneName == "Hnrnpul1" & shRNA == "KD_day12_shSCR")$CPMs 

ggplot(mKDeffic, aes(x=Condition, y=CPMs, color=GeneName)) +
  facet_wrap(.~GeneName)+
  geom_line(size=.4,color="gray80") +
  geom_point() +
  theme_classic() +
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=30,hjust=1,vjust=1),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
ggsave("KDefficiency.pdf", width=7,height=6,device = cairo_pdf)
write.table(mKDeffic,"KDefficiency.txt",sep="\t",row.names = F,quote=F)

## CHECK OE EFFICIENCY
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/GEx_EdgeR/OE/")
OEeffic <- subset(OEs_edgeR, GeneName == "Tia1")       
mOEeffic <- melt(OEeffic, variable.name = "Sample", value.name = "CPMs")
mOEeffic$Condition <- gsub("_1","",mOEeffic$Sample)
mOEeffic$Condition <- gsub("_2","",mOEeffic$Condition)
mOEeffic$Condition <- factor(mOEeffic$Condition, levels =c("OE_day00_Ctrl","OE_day12_Empty","OE_day12_T7TIA1"))
mOEeffic$Replicate <- str_sub(mOEeffic$Sample,-1)
head(mOEeffic)

mOEeffic_av <- aggregate(CPMs ~ ENSGID+GeneName+Condition, mOEeffic, mean)
subset(mOEeffic_av, GeneName == "Tia1" & Condition == "OE_day12_T7TIA1")$CPMs / subset(mOEeffic_av, GeneName == "Tia1" & Condition == "OE_day12_Empty")$CPMs 

ggplot(mOEeffic, aes(x=Condition, y=CPMs, color=GeneName)) +
  facet_wrap(.~GeneName)+
  geom_line(size=.4,color="gray80") +
  geom_point() +
  theme_classic() +
  theme (text = element_text(color="grey20",size=11),
         axis.title = element_text(face="bold"),
         axis.text.x = element_text(angle=30,hjust=1,vjust=1),
         legend.position = "bottom", legend.title = element_text(face="bold"),
         panel.grid.minor.x=element_blank(),
         panel.grid.major.y=element_line(linetype = "dashed",color="gray80"),
         strip.background = element_blank(),strip.text = element_text(face="bold"))
ggsave("OEefficiency.pdf", width=7,height=6,device = cairo_pdf)
# write.table(mOEeffic,"OEefficiency.txt",sep="\t",row.names = F,quote=F)
