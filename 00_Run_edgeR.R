require(edgeR)
require(tidyverse)
require(org.Mm.eg.db)


## FUNCTIONS 
DiffGEx <- function(DGEList, conditions_to_compare,
                    top = 500, diff_correction = "BH", pval = 0.001,
                    write_output = c("yes","no")) {
  et <- exactTest(DGEList, pair = conditions_to_compare)  ## Pair = group to compare
  tt <- topTags(et, top)
  de <- decideTestsDGE(et, adjust.method=diff_correction, p.value=pval)
  print(summary(de))  #-1, 0 and 1 are for down-regulated, non-differentially expressed and up-regulated
  
  tt$table$Group <- rep(paste0("top",top),nrow(tt$table))
  tmp <- merge(x=et$table,y=et$genes,by="row.names")
  tmp2 <- merge(x=tmp,y=tt$table, by=c("Symbol","logFC","logCPM","PValue"),all=T)
  out <- merge(x=tmp2,y=de,by.x="Row.names",by.y="row.names")
  rownames(out) <- out$Row.names
  colnames(out) <- c("ENSGID","Symbol","logFC","logCPM","PValue","FDR","Group","DE")
  
  p <- unlist(str_split(as.character(pval),"[.]"))[2]
  filename <- paste0("DiffExp_",conditions_to_compare[1],"_",conditions_to_compare[2],"_top",top,"_DEpval",p,".txt")
  if(write_output == "yes") {
    write.table(x= out, file = filename,row.names = F, sep="\t",quote =F)} else {return(out)}
}



## READ COUNT FILES FROM STAR
setwd("~/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/GEx_EdgeR/")
files <- list.files (path= "~/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/MAPPING/", pattern = "*ReadsPerGene.out.tab")
files_KD <- files[grep("^KD_",files)]
files_OE <- files[grep("^OE_",files)]
    dataset_KD <- sapply(files_KD, function(x) unlist(strsplit(x,split="_"))[1])
    dataset_OE <- sapply(files_OE, function(x) unlist(strsplit(x,split="_"))[1])
    day_KD <- sapply(files_KD, function(x) unlist(strsplit(x,split="_"))[2])
    day_OE <- sapply(files_OE, function(x) unlist(strsplit(x,split="_"))[2])
    condition_KD <- sapply(files_KD, function(x) unlist(strsplit(x,split="_"))[3])
    condition_OE <- sapply(files_OE, function(x) unlist(strsplit(x,split="_"))[3])
    replicates_KD <- sapply(files_KD, function(x) unlist(strsplit(x,split="_"))[4])
    replicates_OE <- sapply(files_OE, function(x) unlist(strsplit(x,split="_"))[4])
    samples_KD <- paste(day_KD,condition_KD, sep="_")
    samples_OE <- paste(day_OE,condition_OE, sep="_")
labels_KD <- paste(dataset_KD,samples_KD,replicates_KD, sep="_")
labels_OE <- paste(dataset_OE,samples_OE,replicates_OE, sep="_")


## CREATE DGELIST OBJECT FROM MANY FILES CONTAINING COUNTS OF SAMPLES
dgList_KD <- readDGE(skip=4, ## skip first 4 lines cause they contain the summary
                     files = files_KD, path = "/Volumes/jvalcarcel/cvivori/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/MAPPING/", 
                     columns=c(1,2), ## the 2nd column contains the unstranded reads!!!!!!!  
                     group = samples_KD, labels = labels_KD)
dgList_KD
head(dgList_KD$counts)
summary(dgList_KD$counts)
#total number of readcounts per sample       colSums(dgList_KD$counts)

dgList_OE <- readDGE(skip=4, ## skip first 4 lines cause they contain the summary
                     files = files_OE, path = "/Volumes/jvalcarcel/cvivori/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/MAPPING/", 
                     columns=c(1,2), ## the 2nd column contains the unstranded reads!!!!!!!  
                     group = samples_OE, labels = labels_OE)
dgList_OE
head(dgList_OE$counts)
summary(dgList_OE$counts)


## CONVERT ENSGIDs in human-readable GENE NAMES
IDs_KD <- rownames(dgList_KD$counts); IDs_OE <- rownames(dgList_OE$counts)
ENSGIDs_KD <- sapply(IDs_KD, function(x) unlist(strsplit(x,split="[.]"))[1]); ENSGIDs_OE <- sapply(IDs_OE, function(x) unlist(strsplit(x,split="[.]"))[1])
Symbol_KD <- mapIds(org.Mm.eg.db, keys=ENSGIDs_KD, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
Symbol_OE <- mapIds(org.Mm.eg.db, keys=ENSGIDs_OE, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
## ADD IT TO THE LIST
dgList_KD$genes <- data.frame (Symbol = Symbol_KD)
dgList_OE$genes <- data.frame (Symbol = Symbol_OE)
dgList_KD


# SAVE ORIGINAL dgList
# dgList_KD_orig <- dgList_KD
# dgList_KD_orig <- calcNormFactors(dgList_KD_orig)
# dgList_KD_orig$samples
# dgList_OE_orig <- dgList_OE
# dgList_OE_orig <- calcNormFactors(dgList_OE_orig)
# dgList_OE_orig$samples




#### SET FILTERS HERE
# 5 cpms at least in 33% of samples? 5 for KDs and 6 for OE
threshold_cpm <- 5    
threshold_nsamples_KD <- 8    # 33% of samples
threshold_nsamples_OE <- 2    

## Filtering for minimal gene expression
keep_KD <- rowSums(cpm(dgList_KD)>threshold_cpm) >= threshold_nsamples_KD    
keep_OE <- rowSums(cpm(dgList_OE)>threshold_cpm) >= threshold_nsamples_OE    
dgList_KD <- dgList_KD[keep_KD,]
dgList_OE <- dgList_OE[keep_OE,]
## Also, all genes have to have at least 1 cpm
keep_KD <- rowSums(cpm(dgList_KD)>1) >= 12    
keep_OE <- rowSums(cpm(dgList_OE)>1) >= 6   
dgList_KD <- dgList_KD[keep_KD,]
dgList_OE <- dgList_OE[keep_OE,]
dim(dgList_KD_orig); dim(dgList_KD)
dim(dgList_OE_orig); dim(dgList_OE)
## Recalculate the library size accordingly
dgList_KD$samples$lib.size <- colSums(dgList_KD$counts)
dgList_OE$samples$lib.size <- colSums(dgList_OE$counts)


## Calculate Normalization Factor
dgList_KD <- calcNormFactors(dgList_KD); #ps_dgList_KD <- calcNormFactors(ps_dgList_KD)
dgList_KD$samples
dgList_OE <- calcNormFactors(dgList_OE); #ps_dgList_OE <- calcNormFactors(ps_dgList_OE)
dgList_OE$samples


## MD PLOT
setwd("~/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/GEx_EdgeR/KD/")
  filter_title_KD = paste("Filter: cpm >", threshold_cpm, "in", threshold_nsamples_KD, "samples", sep=" ")
  for (x in c(1:ncol(dgList_KD$counts))) {
    sample_name = colnames(dgList_KD$counts)[x]
    pdf(file = paste("MD_cpm",threshold_cpm,"_nsamples",threshold_nsamples_KD,"_",sample_name,".pdf", sep=""), width =10, height=7)
      par(mfrow=c(1,2))
      plotMD(cpm(dgList_KD_orig, log=T),column=x, main = "Before filtering"); abline(h=0,col="red",lty=2,lwd=1)
      plotMD(cpm(dgList_KD, log=T),column=x, main = filter_title_KD ); abline(h=0,col="red",lty=2,lwd=1)
      title(sample_name, outer = TRUE, cex = 1.5, line = -1)
    dev.off()
  }
setwd("~/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/GEx_EdgeR/OE/")
  filter_title_OE = paste("Filter: cpm >", threshold_cpm, "in", threshold_nsamples_OE, "samples", sep=" ")
  for (x in c(1:ncol(dgList_OE$counts))) {
    sample_name = colnames(dgList_OE$counts)[x]
    pdf(file = paste("MD_cpm",threshold_cpm,"_nsamples",threshold_nsamples_OE,"_",sample_name,".pdf", sep=""), width =10, height=7)
    par(mfrow=c(1,2))
    plotMD(cpm(dgList_OE_orig, log=T),column=x, main = "Before filtering"); abline(h=0,col="red",lty=2,lwd=1)
    plotMD(cpm(dgList_OE, log=T),column=x, main = filter_title_OE ); abline(h=0,col="red",lty=2,lwd=1)
    title(sample_name, outer = TRUE, cex = 1.5, line = -1)
    dev.off()
  }

  
### MDS PLOTS
    setwd("~/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/GEx_EdgeR/KD/")
    pdf(file = paste("MDS_cpm",threshold_cpm,"_nsamples",threshold_nsamples_KD,"_KD.pdf", sep=""),width=10,height=7)
    par(mfrow=c(1,2))
      plotMDS(dgList_KD_orig, method="bcv",col=replicates_KD, main = "Before filtering")
      plotMDS(dgList_KD, method="bcv",col=replicates_KD, main = filter_title_KD )
      title("KD", outer = TRUE, cex = 1.5, line = -1)
    dev.off()
    setwd("/Volumes/jvalcarcel/cvivori/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/GEx_EdgeR/OE/")
    pdf(file = paste("MDS_cpm",threshold_cpm,"_nsamples",threshold_nsamples_OE,"_MEFstoiPS.pdf", sep=""),width=10,height=7)
    par(mfrow=c(1,2))
      plotMDS(dgList_OE_orig, method="bcv",col=replicates_OE, main = "Before filtering")
      plotMDS(dgList_OE, method="bcv",col=replicates_OE, main = filter_title_OE )
      title("OE", outer = TRUE, cex = 1.5, line = -1)
    dev.off()

    
## OUTPUT TABLES WITH CPMs 
cpm_KD <- cpm(dgList_KD)
ENSGID_KD <- sapply(rownames(cpm_KD), function(x) unlist(strsplit(x,split="[.]"))[1])
GeneName_KD <- mapIds(org.Mm.eg.db, keys=ENSGID_KD, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
df_KD <- data.frame(ENSGID=ENSGID_KD, GeneName=GeneName_KD, cpm_KD)
head(df_KD)   
# write.table(x= df_KD,file = "FilteredGeneCpm_KD.txt", row.names = F, sep="\t",quote =F)

cpm_OE <- cpm(dgList_OE)
ENSGID_OE <- sapply(rownames(cpm_OE), function(x) unlist(strsplit(x,split="[.]"))[1])
GeneName_OE <- mapIds(org.Mm.eg.db, keys=ENSGID_OE, keytype="ENSEMBL",column="SYMBOL", multiVals = "first")
df_OE <- data.frame(ENSGID=ENSGID_OE,GeneName= GeneName_OE, cpm_OE)
head(df_OE)   
# write.table(x= df_OE,file = "FilteredGeneCpm_OE.txt", row.names = F, sep="\t",quote =F)

    
### ESTIMATE DISPERSION ACCORDING TO EDGER
d1_KD <- estimateCommonDisp(dgList_KD, verbose = T)
d1_OE <- estimateCommonDisp(dgList_OE, verbose = T)
d1_KD <- estimateTagwiseDisp(d1_KD)
d1_OE <- estimateTagwiseDisp(d1_OE)

plotBCV(d1_KD)
plotBCV(d1_OE)

design.mat_KD <- model.matrix(~ 0+ samples_KD) 
rownames(design.mat_KD) <- rownames(dgList_KD$samples)
d2_KD <- estimateGLMCommonDisp(dgList_KD,design.mat_KD)
d2_KD <- estimateGLMTrendedDisp(d2_KD,design.mat_KD, method="auto")
d2_KD <- estimateGLMTagwiseDisp(d2_KD,design.mat_KD)
plotBCV(d2_KD)

design.mat_OE <- model.matrix(~  0+ samples_OE) 
rownames(design.mat_OE) <- rownames(dgList_OE$samples)
d2_OE <- estimateGLMCommonDisp(dgList_OE,design.mat_OE)
d2_OE <- estimateGLMTrendedDisp(d2_OE,design.mat_OE, method="auto")
d2_OE <- estimateGLMTagwiseDisp(d2_OE,design.mat_OE)
plotBCV(d2_OE)




## SET THRESHOLDS TO USE FOR DIFF.GEx
top_value <- 500
pval_value <- 0.001

## CALCULATE DIFFERENTIALLY EXPRESSED GENES in KD
DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day00_Ctrl","day12_shSCR"),
         top = top_value, diff_correction = "BH", pval = pval_value,
         write_output = "yes")

DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day00_Ctrl","day12_shCPSF3n1"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")
DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day00_Ctrl","day12_shCPSF3n5"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")  

DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day00_Ctrl","day12_shUL1n3"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")
DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day00_Ctrl","day12_shUL1n4"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")  
  
DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day12_shSCR","day12_shCPSF3n1"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")
DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day12_shSCR","day12_shCPSF3n5"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")

DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day12_shSCR","day12_shUL1n3"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")
DiffGEx(DGEList = d1_KD, conditions_to_compare = c("day12_shSCR","day12_shUL1n4"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")

## CALCULATE DIFFERENTIALLY EXPRESSED GENES in OE
DiffGEx(DGEList = d1_OE, conditions_to_compare = c("day00_Ctrl","day12_Empty"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")
DiffGEx(DGEList = d1_OE, conditions_to_compare = c("day00_Ctrl","day12_T7TIA1"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")
DiffGEx(DGEList = d1_OE, conditions_to_compare = c("day12_Empty","day12_T7TIA1"),
        top = top_value, diff_correction = "BH", pval = pval_value,
        write_output = "yes")


