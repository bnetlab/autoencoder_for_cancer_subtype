library(CancerSubtypes)
library(RTCGA.mRNA)
rm(list = ls())

data("BRCA.mRNA")
mRNA=t(as.matrix(BRCA.mRNA[,-1]))
colnames(mRNA)=BRCA.mRNA[,1]
data.checkDistribution(mRNA)

# data cleaning
index=which(is.na(mRNA))
result=data.imputation(mRNA,fun="median")
result=data.normalization(mRNA,type="feature_zscore",log2=FALSE)


library("CancerSubtypes")
load("GBM_GeneEXp.rda")
load("GBM_clinical.rda")

# clustering
GeneExp=list(GBM_GeneEXp)
result=ExecuteCC(clusterNum=3,d=GeneExp,maxK=10,clusterAlg="hc",distance="pearson",title="GBM")

#plot silhouette score

sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
plot(sil)

#plot servival analysis
data(time)
data(status)

data1=FSbyCox(GeneExp,time,status,cutoff=0.05)
p_value=survAnalysis(mainTitle="GBM1",time,status,result$group,
                     result$distanceMatrix,similarity=TRUE)

# statistical significance of clustering

sigclust=sigclustTest(GeneExp,result$group, nsim=1000, nrep=1, icovest=1)
sigclust