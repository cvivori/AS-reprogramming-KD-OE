require(tidyverse)
require(gdata)
require(biomaRt)
hmart = useMart (biomart = "ENSEMBL_MART_ENSEMBL",host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl")
require(org.Mm.eg.db)   #keytypes(org.Mm.eg.db)
require(org.Hs.eg.db)   #keytypes(org.Hs.eg.db)


## FUNCTIONs
Human2Mouse <- function (hENSGId_vector)         
{                                    
  IDs2genes = getBM (attributes = c("mmusculus_homolog_ensembl_gene", "ensembl_gene_id"),  # add the attributes you need
                     filters = "ensembl_gene_id",     #type of IDs you want to convert
                     values = hENSGId_vector,            #vector containing your IDs 
                     mart = hmart)
  mENSGId_vector= IDs2genes[,"mmusculus_homolog_ensembl_gene"]
  return(mENSGId_vector)
}


setwd("~/Dropbox (CRG ADV)/Personal_Claudia/Cl@udia/PhD/Data/1912 RNAseq shCPSF3 shUL1 T7TIA1/VAST-TOOLS_2.2.2/OE/GO_analysis/")

## IMPORT REVIGO TABLE - SELECT SECRETION-RELATED CATEGORIES
GO_TIA1dep_REVIGO <- read.csv(file= "REVIGO_GOPROCESS_TIA1dep_vsFilINCL.csv",header=T,fill=NA)
GO_TIA1dep_REVIGO_secretion <- GO_TIA1dep_REVIGO %>% 
  filter(representative == 90277)
secretion_IDs <- GO_TIA1dep_REVIGO_secretion$term_ID

## IMPORT GORILLA TABLE - SELECT SECRETION-RELATED GENES
GO_TIA1dep_PROCESS <- read.xls(xls= "GOPROCESS_TIA1dep_vsFilINCL.xlsx",sheet="GOPROCESS_TIA1dep_vsFilINCL",header=T,fill=NA)
GO_TIA1dep_PROCESS_secretion <- GO_TIA1dep_PROCESS %>% 
  filter(GO.Term %in% secretion_IDs)

## EXTRACT GENES FROM SECRETION SETS
GO_Genes <- unlist(apply(GO_TIA1dep_PROCESS_secretion, 1, function(x) str_split(x["Genes"], pattern = ",")))
GO_Genes <- GO_Genes[grepl(GO_Genes,pattern = " - ")]
GO_Genes <- str_trim(str_split_fixed(GO_Genes, pattern = " - ",n=2)[,1], side = "both")
GO_Genes <- str_replace(GO_Genes, pattern = "\\[",replacement = "")
GO_Genes

## CONVERT TO ENSGIDs
GO_ENSGIDs <- mapIds(org.Mm.eg.db, keys=GO_Genes, column='ENSEMBL', keytype='ALIAS',multiVals="first")	
GO_ENSGIDs <- unique(GO_ENSGIDs)
length(GO_ENSGIDs)

## CONVERT HUMAN SENESCENCE/SASP-RELATED GENES TO MOUSE ENSGIDs
SEN_genes_hs <-as.character(read.csv(file = "GeneCards_SASP_OR_Senescence.csv")$Gene.Symbol)
SEN_ENSGIDs_hs <- mapIds(org.Hs.eg.db, keys=SEN_genes_hs, column='ENSEMBL', keytype='ALIAS',multiVals="first")	
SEN_ENSGIDs_mm <- Human2Mouse(SEN_ENSGIDs_hs) 
SEN_ENSGIDs_mm <- unique(SEN_ENSGIDs_mm[grepl("ENSMUSG",SEN_ENSGIDs_mm)])
length(SEN_ENSGIDs_mm)

## PROPORTION OF CORRESPONDING ENSGIDs
sum(is.element(GO_ENSGIDs,SEN_ENSGIDs_mm))/length(GO_ENSGIDs)

