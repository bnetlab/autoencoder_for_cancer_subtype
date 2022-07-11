# Download miRNA and DNA mathelation data

rm(list=ls())
library(CancerSubtypes)
library("TCGAbiolinks")
library("SummarizedExperiment")

cancerType <- "GBM"
directory <- "./GDC/"
CancerProject <- paste0("TCGA-",cancerType)
DataDirectory <- paste0(directory,"GDC_",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","DNA",".rda")

## DNA data
query <- GDCquery(project = CancerProject,
                  data.category = "DNA methylation", 
                  platform = "Illumina Human Methylation 27", 
                  legacy = TRUE)
query_case1 <- query$results[[1]]$cases

queryDown<- GDCquery(project = CancerProject,
                  data.category = "DNA methylation", 
                  platform = "Illumina Human Methylation 27", 
                  barcode = query_case1,
                  legacy = TRUE)

GDCdownload(queryDown, directory = DataDirectory)
data <- GDCprepare(query = queryDown,
                    save = TRUE,
                    directory =  DataDirectory,
                    save.filename = paste0(DataDirectory, "_","Methylation_27",".rda"),
                    summarizedExperiment = TRUE)
data2 <- assay(data, 1)

