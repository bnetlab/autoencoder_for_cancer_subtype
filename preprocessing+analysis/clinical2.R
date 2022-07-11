rm(list=ls())
library(CancerSubtypes)
library("TCGAbiolinks")

#loading expresion data
geneExpresion <- read.csv('../data/TCGA/TCGA-BRCA-tumor.csv', header=TRUE, sep=',', row.names = 'index')

# Downloading clinical data for a cancer type
cancerType <- "BRCA"
directory <- "./BRCA/"
CancerProject <- paste0("TCGA-",cancerType)
DataDirectory <- paste0(directory,"BRCA_",gsub("-","_",CancerProject))
query_clin <- GDCquery(project = CancerProject,file.type = "xml",data.category = "Clinical")
clinical_case <- query_clin$results[[1]]$cases
queryDown_clin <- GDCquery(project = CancerProject,data.category = "Clinical",barcode = clinical_case ,file.type = "xml")
GDCdownload(queryDown_clin)

# clinical data cleaning
clinical <- GDCprepare_clinic(queryDown_clin,clinical.info = "patient")
clinical <- clinical[!duplicated(clinical$bcr_patient_barcode),]
rownames(clinical) <- clinical$bcr_patient_barcode
GBM_clinical=clinical[,c("days_to_death","days_to_last_followup","vital_status")]
index4=which(is.na(GBM_clinical[,"days_to_death"]))
GBM_clinical[index4,"days_to_death"]=GBM_clinical[index4,"days_to_last_followup"]
status=as.vector(GBM_clinical[,"vital_status"])
status[status=="Alive"]=0
status[status=="Dead"]=1
GBM_clinical=cbind(GBM_clinical,"status"=as.numeric(status))
colnames(GBM_clinical)[1]="time"

# intersect sample
intersect_sample <- intersect(substr(rownames(geneExpresion),1,12), rownames(clinical))
index5=match(intersect_sample,substr(rownames(geneExpresion),1,12))
index6=match(intersect_sample,rownames(GBM_clinical))
geneExpresion_data_main=geneExpresion[index5,]
GBM_clinical=GBM_clinical[index6,]

# Clustering using spactral clustering
geneExpresion_data = t(log2(geneExpresion_data_main+1))
geneExpresion_data = data.normalization(geneExpresion_data ,type="feature_zscore",log2=FALSE)
GBM <- list(GeneExp=geneExpresion_data)
result = ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
group <- result$group
p_value=survAnalysis(mainTitle="GBMOriginal",GBM_clinical$time,GBM_clinical$status,group,
                     distanceMatrix=result$distanceMatrix,similarity=TRUE)

