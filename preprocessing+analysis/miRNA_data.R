# Download miRNA and DNA mathelation data

rm(list=ls())
library(CancerSubtypes)
library("TCGAbiolinks")
library("SummarizedExperiment")

cancerType <- "GBM"
directory <- "./GDC/"
CancerProject <- paste0("TCGA-",cancerType)
DataDirectory <- paste0(directory,"GDC_",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","miRNA",".rda")

query1 <- GDCquery(project = CancerProject,
                   data.category = "Gene expression",
                   data.type = "miRNA gene quantification",
                   legacy = TRUE)
query_case1 <- query1$results[[1]]$cases

queryDown1 <- GDCquery(project = CancerProject,
                       data.category = "Gene expression",
                       data.type = "miRNA gene quantification",
                       barcode = query_case1,
                       legacy = TRUE)

GDCdownload(queryDown1,directory = DataDirectory)

dataPrep1 <- GDCprepare(query = queryDown1,
                        save = TRUE,
                        directory =  DataDirectory,
                        save.filename = paste0(DataDirectory, "_","miRNA",".rda"),
                        summarizedExperiment = FALSE)
