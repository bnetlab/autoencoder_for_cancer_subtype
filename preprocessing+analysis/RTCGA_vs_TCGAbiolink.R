rm(list=ls())

# mRNA data from RTCGA
library("RTCGA.mRNA")
data(BRCA.mRNA)
mRNA=t(as.matrix(BRCA.mRNA[,-1]))
colnames(mRNA)=BRCA.mRNA[,1]
dim(mRNA) 

index1=which(as.numeric(substr(colnames(mRNA),14,15))>9)
mRNA_Normal=mRNA[,index1]
dim(mRNA_Normal) 
index2=which(as.numeric(substr(colnames(mRNA),14,15))<9)
mRNA_Case=mRNA[,index2]
dim(mRNA_Case) 

# mRNA data from TCGA biolink

library("TCGAbiolinks")
library("SummarizedExperiment")

cancerType <- "BRCA"
directory <- "./BRCA/"
CancerProject <- paste0("TCGA-",cancerType)
DataDirectory <- paste0(directory,"GDC_",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","AgilentG4502A_07_1",".rda")

######GBM Gene expression Data1: AgilentG4502A_07_1########

query1 <- GDCquery(project = CancerProject,
                   data.category = "Gene expression",
                   data.type = "Gene expression quantification",
                   platform = "Illumina HiSeq",
                   legacy = TRUE)
query_case1 <- query1$results[[1]]$cases

queryDown1 <- GDCquery(project = CancerProject,
                       data.category = "Gene expression",
                       data.type = "Gene expression quantification",
                       platform = "Illumina HiSeq",
                       barcode = query_case1,
                       legacy = TRUE)

GDCdownload(queryDown1,directory = DataDirectory)

dataPrep1 <- GDCprepare(query = queryDown1,
                        save = TRUE,
                        directory =  DataDirectory,
                        save.filename = paste0(DataDirectory, "_","Illumina HiSeq",".rda"),
                        summarizedExperiment = TRUE)
data1 <- assay(dataPrep1, 1)
##data imputation for missing measurements
data1=data.imputation(data1, fun = "microarray")

######GBM Gene expression Data2: AgilentG4502A_07_2########

query2 <- GDCquery(project = CancerProject,
                   data.category = "Gene expression",
                   data.type = "Gene expression quantification",
                   platform = "AgilentG4502A_07_3",
                   legacy = TRUE)
query_case2 <- query2$results[[1]]$cases

queryDown2 <- GDCquery(project = CancerProject,
                       data.category = "Gene expression",
                       data.type = "Gene expression quantification",
                       platform = "AgilentG4502A_07_3",
                       barcode = query_case2,
                       legacy = TRUE)
GDCdownload(queryDown2,directory = DataDirectory)

dataPrep2 <- GDCprepare(query = queryDown2,
                        save = TRUE,
                        directory =  DataDirectory,
                        save.filename = paste0(DataDirectory, "_","AgilentG4502A_07_3",".rda"))

data2 <- assay(dataPrep2, 1)
data2=data.imputation(data2, fun = "microarray")

##combined two platform
GBM_mRNA=cbind(data1,data2)

###Extract the normal samples
index1=which(as.numeric(substr(colnames(GBM_mRNA),14,15))>9)
GBM_mRNA_Normal=GBM_mRNA[,index1]

###Extract the PRIMARY SOLID TUMOR("TP") samples
index2=which(substr(colnames(GBM_mRNA),14,15)=="01")
GBM_mRNA_Tumor=GBM_mRNA[,index2]

## observation:
# BRCA RTCGA = AgilentG4502A_07_3

##################miRNA#######################
# # miRNA data from RTCGA
# library("RTCGA.miRNASeq")
# miRNA=t(as.matrix(COAD.miRNASeq[,-1]))
# tokeep <- seq(1, ncol(miRNA), 3)
# miRNA <-  miRNA[,tokeep ]
# miRNA <- miRNA[-1,]
# #result3=data.normalization(res3,type="feature_zscore",log2=TRUE)
# ###Extract the normal samples
# index1=which(as.numeric(substr(colnames(miRNA),14,15))>9)
# miRNA_Normal=miRNA[,index1]
# ###Extract the case samples
# index2=which(as.numeric(substr(colnames(miRNA),14,15))<9)
# miRNA_Case=miRNA[,index2]
# dim(miRNA)
# dim(miRNA_Normal)
# dim(miRNA_Case)
# 
# #miRNA data TCGAbiolink
# 
# library("TCGAbiolinks")
# library("SummarizedExperiment")
# 
# cancerType <- "BRCA"
# directory <- "./BRCA/"
# CancerProject <- paste0("TCGA-",cancerType)
# DataDirectory <- paste0(directory,"GDC_",gsub("-","_",CancerProject))
# FileNameData <- paste0(DataDirectory, "_","miRNA",".rda")
# 
# query1 <- GDCquery(project = CancerProject,
#                    data.category = "Gene expression",
#                    data.type = "miRNA gene quantification",
#                    legacy = TRUE)
# query_case1 <- query1$results[[1]]$cases
# 
# queryDown1 <- GDCquery(project = CancerProject,
#                        data.category = "Gene expression",
#                        data.type = "miRNA gene quantification",
#                        barcode = query_case1,
#                        legacy = TRUE)
# 
# GDCdownload(queryDown1,directory = DataDirectory)
# 
# dataPrep1 <- GDCprepare(query = queryDown1,
#                         save = TRUE,
#                         directory =  DataDirectory,
#                         save.filename = paste0(DataDirectory, "_","miRNA",".rda"),
#                         summarizedExperiment = FALSE)

