library(data.table)
library(dplyr)
library(trqwe)

source("./counts_to_TPM.R")
source("./Imputation.R")
source("./EMT_score_func.R")

## Folder with all the raw read counts
setwd("./Data")

fileList = list.files(pattern = "*_counts.rds")
fileList = sapply(strsplit(fileList, split = "_"), function(x) x[1])

corMat = matrix(0, nrow = length(fileList), ncol = 3)
colnames(corMat) = c("Sample","76GS-KS_Cor","76GS-KS_Pval")

for(dataNum in 1:length(fileList)){
    countToTpm_sc(fileList[dataNum])
    impute(fileList[dataNum])
    counts = mcreadRDS(paste0("../Data_generated/",fileList[dataNum],"_imputed_exact.rds"), mc.cores=4)
    counts = rnaToMA(counts)
    EMT76GS_Score = EMT76GS(counts)
    KS_score=KSScore(counts)
    MLR_score = mlrEMTPred(counts,fileList[dataNum])
    writeEMTscore(fileList[dataNum], EMT76GS_Score, KS_score)
    corMat[dataNum, ] = c(fileList[dataNum],all_scoreCor(list(EMT76GS_Score[,1], KS_score[,1])))

}

fwrite(corMat, "../Output/EMT_Score_Correlations.tsv",row.names=F,sep='\t')
