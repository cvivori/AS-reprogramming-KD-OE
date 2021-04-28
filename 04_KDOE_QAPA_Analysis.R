require(tidyverse)
require(ggrepel)
require(reshape2)
require(viridis)
require(VennDiagram)
require(org.Mm.eg.db)
require(bedr)


## IMPORT PAU TABLE FROM QAPA
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/QAPA/")
PAU <- read.table("pau_results.txt",header=T,sep="\t",fill=NA)
rownames(PAU)=PAU$APA_ID
head(PAU)
dim(PAU)

# Remove transcripts with only one polyA and expressed 3TPM in less than 83% of samples
tofilter_PAU <- subset(PAU, Num_Events >1); dim(tofilter_PAU)
tofilter_TPM <- tofilter_PAU[,grepl(".TPM",colnames(tofilter_PAU))]
tofilter_KD_TPM <- tofilter_TPM[,grepl("KD_", colnames(tofilter_TPM))]
tofilter_OE_TPM <- tofilter_TPM[,grepl("OE_", colnames(tofilter_TPM))]
threshold <- c(KD = 10, OE = 5) ## 83%
filtered_PAU <- list(KD = tofilter_PAU[which(rowSums(tofilter_KD_TPM>3)>threshold["KD"]),!grepl(".TPM", colnames(tofilter_PAU))],
                     OE = tofilter_PAU[which(rowSums(tofilter_OE_TPM>3)>threshold["OE"]),!grepl(".TPM", colnames(tofilter_PAU))])
lapply(filtered_PAU, dim)
lapply(filtered_PAU, head)

## Extract samples and controls
samples <- unique(str_remove(colnames(filtered_PAU$KD)[-c(1:12)],pattern = "_[12].PAU"))
KD_samples <- samples[grepl("KD_",samples)]
ctrl_cols <- c("KD_day00_Ctrl_1.PAU","KD_day00_Ctrl_2.PAU")
KD_PAU_av = KD_PAU_range = data.frame(APA_ID = filtered_PAU$KD$APA_ID)
rownames(KD_PAU_av) <- KD_PAU_av$APA_ID; rownames(KD_PAU_range) <- KD_PAU_range$APA_ID

## Calculate average and range
for (sample in KD_samples) {
cols <- colnames(filtered_PAU$KD)[grep(sample,colnames(filtered_PAU$KD))]
  KD_PAU_av[,sample] = apply(filtered_PAU$KD, 1, function(x) sum(as.numeric(x[cols]))/length(cols) )
  KD_PAU_range[,sample] = apply(filtered_PAU$KD, 1, function(x) min(abs( c( as.numeric(x[cols])-as.numeric(x[ctrl_cols]), as.numeric(x[rev(cols)])-as.numeric(x[ctrl_cols])) )) )
}
head(KD_PAU_av); head(KD_PAU_range)

## Generate list containing all values and calculate ∆PAU
all_KD <- colnames(filtered_PAU$KD[grepl("KD_",colnames(filtered_PAU$KD))])
PAU_list <- list (unfiltered_PAU = tofilter_PAU[,c("APA_ID","Gene","Gene_Name",all_KD)],
              filtered_PAU = filtered_PAU$KD[,c("APA_ID","Gene","Gene_Name",all_KD)],
              PAU_av = KD_PAU_av,
              PAU_range = KD_PAU_range,
              dPAU = data.frame(APA_ID = KD_PAU_av$APA_ID,
                      KD_day12_shSCR = KD_PAU_av$KD_day12_shSCR - KD_PAU_av$KD_day00_Ctrl,
                      KD_day12_shCPSF3n1 = KD_PAU_av$KD_day12_shCPSF3n1 - KD_PAU_av$KD_day00_Ctrl,
                      KD_day12_shCPSF3n5 = KD_PAU_av$KD_day12_shCPSF3n5 - KD_PAU_av$KD_day00_Ctrl,
                      KD_day12_shUL1n3 = KD_PAU_av$KD_day12_shUL1n3 - KD_PAU_av$KD_day00_Ctrl,
                      KD_day12_shUL1n4 = KD_PAU_av$KD_day12_shUL1n4 - KD_PAU_av$KD_day00_Ctrl ))
## Calculate ∆∆PAU
PAU_list$ddPAU = data.frame(APA_ID = PAU_list$dPAU$APA_ID,
                            KD_day12_shCPSF3n1 = PAU_list$dPAU$KD_day12_shCPSF3n1 - PAU_list$dPAU$KD_day12_shSCR,
                            KD_day12_shCPSF3n5 = PAU_list$dPAU$KD_day12_shCPSF3n5 - PAU_list$dPAU$KD_day12_shSCR,
                            KD_day12_shUL1n3 = PAU_list$dPAU$KD_day12_shUL1n3 - PAU_list$dPAU$KD_day12_shSCR,
                            KD_day12_shUL1n4 = PAU_list$dPAU$KD_day12_shUL1n4 - PAU_list$dPAU$KD_day12_shSCR) 
rownames(PAU_list$dPAU) <- rownames(KD_PAU_av)
rownames(PAU_list$ddPAU) <- rownames(KD_PAU_av)
summary(PAU_list)
head(PAU_list$PAU_av)
head(PAU_list$ddPAU)

## Differential PAU (dPAU >= 20, range >= 5)
dPAU_20 <- list(shSCR = PAU_list$dPAU[which(abs(PAU_list$dPAU$KD_day12_shSCR)>=20),"APA_ID"],
                shC1 = PAU_list$dPAU[which(abs(PAU_list$dPAU$KD_day12_shCPSF3n1)>=20),"APA_ID"],
                shC5 = PAU_list$dPAU[which(abs(PAU_list$dPAU$KD_day12_shCPSF3n5)>=20),"APA_ID"],
                shU3 = PAU_list$dPAU[which(abs(PAU_list$dPAU$KD_day12_shUL1n3)>=20),"APA_ID"],
                shU4 = PAU_list$dPAU[which(abs(PAU_list$dPAU$KD_day12_shUL1n4)>=20),"APA_ID"])
dPAU_20_range5 <- list(shSCR = intersect(dPAU_20$shSCR, PAU_list$PAU_range[which(abs(PAU_list$PAU_range$KD_day12_shSCR)>=5),"APA_ID"]),
                       shC1 = intersect(dPAU_20$shC1, PAU_list$PAU_range[which(abs(PAU_list$PAU_range$KD_day12_shCPSF3n1)>=5),"APA_ID"]),
                       shC5 = intersect(dPAU_20$shC5, PAU_list$PAU_range[which(abs(PAU_list$PAU_range$KD_day12_shCPSF3n5)>=5),"APA_ID"]),
                       shU3 = intersect(dPAU_20$shU3, PAU_list$PAU_range[which(abs(PAU_list$PAU_range$KD_day12_shUL1n3)>=5),"APA_ID"]),
                       shU4 = intersect(dPAU_20$shU4, PAU_list$PAU_range[which(abs(PAU_list$PAU_range$KD_day12_shUL1n4)>=5),"APA_ID"]) )

# Intersection between shRNAs
dPAU_20_range5$shCPSF3 <- intersect(dPAU_20_range5$shC1,dPAU_20_range5$shC5)
dPAU_20_range5$shUL1 <- intersect(dPAU_20_range5$shU3,dPAU_20_range5$shU4)
## ~40 overlap between shRNAs
length(dPAU_20_range5$shCPSF3) / length(dPAU_20_range5$shC5)
length(dPAU_20_range5$shUL1) / length(dPAU_20_range5$shU4)
# Number of events
unlist(lapply(dPAU_20,length))
unlist(lapply(dPAU_20_range5,length))
# Number of genes
unlist(lapply(dPAU_20_range5, function(x) length(unique(filtered_PAU$KD[which(filtered_PAU$KD$APA_ID %in% x),"Gene"]))))

## CPSF3 and UL1 -dependent polyA events
ddPAU_10 <- list( shC1 = PAU_list$ddPAU[which(abs(PAU_list$ddPAU$KD_day12_shCPSF3n1)>=10),"APA_ID"],
                  shC5 = PAU_list$ddPAU[which(abs(PAU_list$ddPAU$KD_day12_shCPSF3n5)>=10),"APA_ID"],
                  shU3 = PAU_list$ddPAU[which(abs(PAU_list$ddPAU$KD_day12_shUL1n3)>=10),"APA_ID"],
                  shU4 = PAU_list$ddPAU[which(abs(PAU_list$ddPAU$KD_day12_shUL1n4)>=10),"APA_ID"])
# Intersection between shRNAs
ddPAU_10$shCPSF3 <- intersect(ddPAU_10$shC1,ddPAU_10$shC5)
ddPAU_10$shUL1 <- intersect(ddPAU_10$shU3,ddPAU_10$shU4)
## ~44/55overlap between shRNAs
length(ddPAU_10$shCPSF3) / length(ddPAU_10$shC5)
length(ddPAU_10$shUL1) / length(ddPAU_10$shU4)
# Number of events
unlist(lapply(ddPAU_10,length))
# Number of genes
unlist(lapply(ddPAU_10, function(x) length(unique(filtered_PAU$KD[which(filtered_PAU$KD$APA_ID %in% x),"Gene"]))))




## Overlap among APA and AS genes 
AS_genes <- list(CPSF3_dep = Final_KD$CPSF3_dep$GENE,
                 UL1_dep = Final_KD$UL1_dep$GENE)
APA_genes <- list(CPSF3_dep = mapIds(org.Mm.eg.db, 
                                     keys= as.vector(PAU_list$filtered_PAU[which(PAU_list$filtered_PAU$APA_ID %in% ddPAU_10$shCPSF3),"Gene"]), 
                                     keytype="ENSEMBL",column="SYMBOL", multiVals = "first"),
                  UL1_dep = mapIds(org.Mm.eg.db, 
                                     keys= as.vector(PAU_list$filtered_PAU[which(PAU_list$filtered_PAU$APA_ID %in% ddPAU_10$shUL1),"Gene"]), 
                                     keytype="ENSEMBL",column="SYMBOL", multiVals = "first"))
lapply(AS_genes, length)
lapply(APA_genes, length)

length(intersect(AS_genes$CPSF3_dep, APA_genes$CPSF3_dep)); length(intersect(AS_genes$CPSF3_dep, APA_genes$CPSF3_dep))/length(APA_genes$CPSF3_dep)
length(intersect(AS_genes$UL1_dep, APA_genes$UL1_dep)); length(intersect(AS_genes$UL1_dep, APA_genes$UL1_dep))/length(APA_genes$UL1_dep)


## DRAW VENN DIAGRAMs
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/QAPA/")
ovl_genes <- list(CPSF3_dep = calculate.overlap(list(AS_genes$CPSF3_dep, APA_genes$CPSF3_dep)),
                  UL1_dep = calculate.overlap(list(AS_genes$UL1_dep, APA_genes$UL1_dep)))
for (t in names(ovl_genes)) {
  pdf(file=paste("VennDiagram_",t,"_AS_APA_genes.pdf",sep = ""),width = 8, height=8)
draw.pairwise.venn( area1 = length(ovl_genes[[t]]$a1),
                    area2 = length(ovl_genes[[t]]$a2),
                    cross.area = length(ovl_genes[[t]]$a3),
                    category = c("AS genes","APA genes"),
                    lwd=0,
                    fill=c("darkslategray","maroon"),
                    cat.fontface = "bold",
                    cat.fontfamily = "Helvetica",
                    cat.pos =0,
                    #print.mode=c('raw','percent'),
                    fontfamily = "Helvetica"
)
dev.off()
}

## Overlap of APA and AS events 
AS_coords <- list(CPSF3_dep = Final_KD$CPSF3_dep,
                 UL1_dep = Final_KD$UL1_dep)

APA_coords <- list(CPSF3_dep = PAU[which(PAU$APA_ID %in% ddPAU_10$shCPSF3),c("Chr","LastExon.Start","LastExon.End", "APA_ID","Num_Events","Strand")],
                  UL1_dep = PAU[which(PAU$APA_ID %in% ddPAU_10$shUL1),c("Chr","LastExon.Start","LastExon.End", "APA_ID","Num_Events","Strand")])

lapply(AS_coords, dim)
lapply(APA_coords, dim)


## Export to run Bedtools intersect
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/QAPA/")
for (n in names(AS_coords)) {
    write.table(AS_coords[[n]], file=paste("AS_",n,".txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
}
for (n in names(APA_coords)) {
  write.table(APA_coords[[n]], file=paste("APA_",n,".txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = F)
}


## OVERLAP WITH BEDTOOLS
# 4 AS_APA_CPSF3_overlap.txt
# 2 AS_APA_UL1_overlap.txt

#DRAW VENN DIAGRAM
setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/QAPA/")
ovl_coords <- list(CPSF3_dep = 4,
                  UL1_dep = 2)
for (t in names(ovl_coords)) {
  pdf(file=paste("VennDiagram_",t,"_AS_APA_events.pdf",sep = ""),width = 8, height=8)
  draw.pairwise.venn( area1 = nrow(AS_coords[[t]]),
                      area2 = nrow(APA_coords[[t]]),
                      cross.area = ovl_coords[[t]],
                      category = c("AS events","APA events"),
                      lwd=0,
                      fill=c("darkslategray","maroon"),
                      cat.fontface = "bold",
                      cat.fontfamily = "Helvetica",
                      cat.pos =0,
                      #print.mode=c('raw','percent'),
                      fontfamily = "Helvetica"
  )
  dev.off()
}

