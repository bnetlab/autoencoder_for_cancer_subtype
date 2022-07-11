rm(list=ls())
library(CancerSubtypes)
library("RTCGA.mRNA")

# read data and perform DE
res1 <- read.csv(file ="/home/ranap/benchmarkingAMDClassification/data/TCGA/TCGA-BRCA-tumor.csv")
res1_Disease=t(as.matrix(res1[,-1]))
colnames(res1_Disease)=res1[,1]
res1_Disease <- log2(res1_Disease+1)

res2 <- read.csv(file ="/home/ranap/benchmarkingAMDClassification/data/TCGA/TCGA-BRCA-normal.csv")
res1_Normal=t(as.matrix(res2[,-1]))
colnames(res1_Normal)=res2[,1]
res1_Normal <- log2(res1_Normal+1)

res=DiffExp.limma(Tumor_Data=res1_Disease,Normal_Data=res1_Normal,group=NULL,topk=NULL,RNAseq=FALSE)
DE=res[[1]]


## subtype detection

# subtype DTA
rm(list=ls())
library(CancerSubtypes)
library(SNFtool)

cancerType <- "BRCA_Multiview" #change it

mRNA <- readRDS(file="BRCA_mRNA_data.rds")
DNA <- readRDS(file="BRCA_DNA_data.rds")
miRNA <- readRDS(file="BRCA_miRNA_data.rds")
clinical <- readRDS(file="BRCA_clinical_data.rds")
# prepocess, creating PEEP and run DTA
mRNA_Z <- t(scale(t(mRNA)))
#mRNA_Z <- mRNA
mRNA_P <- matrix(0, nrow(mRNA_Z), ncol(mRNA_Z) )
mRNA_P[mRNA_Z>2.5 | mRNA_Z< -2.5]=1

# For DTA feature selection, skip otherwise
write.table(t(mRNA_P), file = "../DTA-Code/DTAin.csv", row.names=FALSE, col.names = FALSE)
system("cd ../DTA-Code && ./run.sh; cd -")
mRNA_feature20 <- scan('../DTA-Code/DTAin_15')
mRNA_feature20
# 2269 10625 12464 12709 13103 13306 13823 14929 16061 16691
# 2269  3798 10625 12464 13103 13306 13823 13960 16061 16691
#  2269  3798 10625 12464 13103 13306 13960 14929 16061 16691


DNA_Z <- t(scale(t(DNA)))
DNA_P <- matrix(0, nrow(DNA_Z), ncol(DNA_Z) )
DNA_P[DNA_Z>2.5 | DNA_Z< -2.5]=1
write.table(t(DNA_P), file = "../DTA-Code/DTAin.csv", row.names=FALSE, col.names = FALSE)
system("cd ../DTA-Code && ./run.sh; cd -")
DNA_feature20 <- scan('../DTA-Code/DTAin_15')
DNA_feature20
# confirmed consistent feature COAD

miRNA_Z <- t(scale(t(miRNA)))
miRNA_P <- matrix(0, nrow(miRNA_Z), ncol(miRNA_Z) )
miRNA_P[miRNA_Z>2.5 | miRNA_Z< -2.5]=1
write.table(t(miRNA_P), file = "../DTA-Code/DTAin.csv", row.names=FALSE, col.names = FALSE)
system("cd ../DTA-Code && ./run.sh; cd -")
miRNA_feature20 <- scan('../DTA-Code/DTAin_15')
miRNA_feature20
# confirmed consistent feature COAD

data1 <- mRNA[mRNA_feature20+1,]
data2 <- DNA[DNA_feature20+1,]
data3 <- miRNA[miRNA_feature20+1,]

cluster=4
data  <- list(GeneExp=data1, DNAmethy=data2, miRNAExp=data3)
result =ExecuteSNF(data, clusterNum=cluster, K=20, alpha=0.5, t=20)
group=result$group
distanceMatrix=result$distanceMatrix
p_value=survAnalysis(mainTitle="Original", clinical$time, clinical$status, group, 
                     distanceMatrix=distanceMatrix, similarity=TRUE)


temp <-cbind(colnames(data1), group)
write.csv(temp, 'BRCA_culster_4_group_k15.csv')
## DE analysis
## reading multiview data
library("RTCGA.mRNA")
data(BRCA.mRNA)
#mRNA data
mRNA_org=t(as.matrix(BRCA.mRNA[,-1])) 
colnames(mRNA_org)=BRCA.mRNA[,1]
mRNA_org=mRNA_org[rowSums(is.na(mRNA_org)) != ncol(mRNA_org), ]
res1=data.imputation(mRNA_org,fun="microarray")
result1=res1
#result1=data.normalization(res1,type="feature_zscore",log2=FALSE)
###Extract the normal samples
index1=which(as.numeric(substr(colnames(result1),14,15))>9)
mRNA_Normal=result1[,index1]

result2=DiffExp.limma(Tumor_Data=mRNA,Normal_Data=mRNA_Normal,group=group,topk=NULL,RNAseq=FALSE)
head(result2[[1]])
head(result2[[2]])
head(result2[[3]])


#Extract top 1500 differentially expressed genes in each subtypes
Subtype1_gene=as.vector(na.omit(result2[[1]]$ID[1:1000]))
Subtype2_gene=as.vector(na.omit(result2[[2]]$ID[1:1000]))
Subtype3_gene=as.vector(na.omit(result2[[3]]$ID[1:1000]))
Subtype4_gene=as.vector(na.omit(result2[[4]]$ID[1:1000]))

library(VennDiagram)
library(VennDiagram)
venn.diagram(filename="BRCA_limma.png",height=300,width=330,list(Subytpe1=Subtype1_gene,Subytpe2=Subtype2_gene,Subytpe3=Subtype3_gene),fill=c("red","green","blue"),alpha=c(0.5,0.5,0.5), cex=1, cat.fontface=4,cat.col = c("dodgerblue", "goldenrod1", "seagreen3"),cat.cex = 1,margin = 0.1,fontfamily=2)

#Functional analysis

# Load libraries
functional <- function(result2) {

library(org.Hs.eg.db)
library(AnnotationDbi)

# Return the Ensembl IDs of gene names
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = result2[[1]]$ID,  # data to use for retrieval
                                           columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retreive for given data
                                           keytype = "SYMBOL") # type of data given in 'keys' argument

# Determine the indices for the non-duplicated genes
# Return only the non-duplicated genes using indices
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

library(DOSE)
library(pathview)
library(clusterProfiler)
library(tidyverse)

## Merge the annotations with the results 
res_ids <- inner_join(result2[[1]], annotations_orgDb, by=c("ID"="SYMBOL")) 

## Create background geneset
## Extract significant geneset
allOE_genes <- as.character(res_ids$ENSEMBL)
sigOE <- dplyr::filter(res_ids, adj.P.Val < 0.01)
sigOE_genes <- as.character(sigOE$ENSEMBL)

## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
dotplot(ego, showCategory=50)
emapplot(ego, showCategory = 50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)
}

functional(result2)
