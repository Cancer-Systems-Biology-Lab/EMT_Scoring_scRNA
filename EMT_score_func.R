library(readxl)
library(matlabr)
library(data.table)

EMT76GS = function(log2Exp){


cat("Calculating EMT Score by using 76 gene signatures ....\n")
   

  remIdx = which(apply(log2Exp,1,function(x) any(x == NaN | x == -Inf)) == TRUE)
  if(length(remIdx) > 0 ) log2Exp = log2Exp[-remIdx, ]

  log2Exp[is.na(log2Exp)] = 0


  sampleNum = ncol(log2Exp)
  genes = rownames(log2Exp)
  exp = apply(log2Exp,2,as.numeric)

  EMTSignature = data.frame(read_excel("../Gene_signatures/76GS/EMT_signature_76GS.xlsx",col_names = TRUE))
  EMTIdx = match(unique(na.omit(EMTSignature[,2])),genes)
  geneFound = length(na.omit(EMTIdx))
  cat(paste(geneFound ,"gene's expression values found \n",sep = " "))
  
  ## Add sd
  sdVal = rnorm(ncol(exp), mean=0, sd= 0.01)
 
  EMTMat = t(apply(exp[na.omit(EMTIdx),], 1, function(x)  x + sdVal))
  row.names(EMTMat) = genes[na.omit(EMTIdx)]


 
  ## get the weights for each gene
  ecadhExp = grep("^CDH1$",row.names(EMTMat))
  if(length(ecadhExp) == 0 ){
    cat("CDH1 gene not found- 76 GS EMT score cannot be calculated\n") 
    EMTScoreStd = rep(0, ncol(exp))
  } else{
    ecadhExpVal = EMTMat[ecadhExp, ]
    weightVal = apply(EMTMat,1,function(x) cor(x,ecadhExpVal))
    EMTMatWt = weightVal*EMTMat
    EMTScore = apply(EMTMatWt,2,sum)
    EMTScoreMean = mean(EMTScore)
    EMTScoreStd = EMTScore-EMTScoreMean

  }  
  emtWrite = cbind(EMTScoreStd)
  colnames(emtWrite) = c("GS76")
  return(emtWrite)

}


KSScore = function(expMat){

  cat("Calculating EMT Score by KS score method ....")

  genes = rownames(expMat)
  exp = apply(expMat,2,as.numeric)
  
  EMTSignature = data.frame(read_excel("../Gene_signatures/KS/EM_gene_signature_cellLine_KS.xlsx",col_names = FALSE))
  commonSig = intersect(EMTSignature[,1],genes)
  EMTExpIdx = match(commonSig,genes)
  EMTExp = exp[EMTExpIdx, ]
  EMTGrpIdx = match(commonSig,EMTSignature[,1])
  geneCat = EMTSignature[EMTGrpIdx,2]
  epiIdx = which(geneCat == "Epi")
  mesIdx = which(geneCat == "Mes")

    ## Perform KS test
    sampleScore2 = matrix(0,nrow=ncol(EMTExp),ncol=6)
    rownames(sampleScore2) = colnames(EMTExp)
    for(i in 1:ncol(EMTExp)){
        ## Two sided test
        ksTwoSided =  ks.test(EMTExp[mesIdx,i],EMTExp[epiIdx,i])
        ## One sided test: ecdf(Mes) > ecdf(Epi)
        ksResGrt = ks.test(EMTExp[mesIdx,i],EMTExp[epiIdx,i],alternative = "greater")
        ## One sided test: ecdf(Epi) > ecdf(Mes)
        ksResLess = ks.test(EMTExp[epiIdx,i],EMTExp[mesIdx,i],alternative = "greater")
        sampleScore2[i, ] = c(ksTwoSided$statistic,ksTwoSided$p.value,
                  ksResGrt$statistic,ksResGrt$p.value,
                  ksResLess$statistic,ksResLess$p.value)
    }

    ## Assign signs to EMT score of sample based on test statistic and pvalue
    finalScore = matrix(0,nrow = nrow(sampleScore2),ncol = 1)
    for (i in 1:nrow(sampleScore2)) {
        if(sampleScore2[i,4] < 0.05){
          finalScore[i, ] = c(-1 * sampleScore2[i,3])
        } else if (sampleScore2[i,6] < 0.05){
          finalScore[i, ] = sampleScore2[i,5]
        } else {

          if(sampleScore2[i,5] == max(c(sampleScore2[i,3],sampleScore2[i,5]))){
            finalScore[i, ] = max(c(sampleScore2[i,3],sampleScore2[i,5]))
          } else {
            finalScore[i, ] = (-1 * max(c(sampleScore2[i,3],sampleScore2[i,5])))
          }
          
      }
  }


  ksOut = cbind(finalScore)
  colnames(ksOut) = c("KS_score")
  rownames(ksOut) = colnames(EMTExp)
  return(ksOut)
}

all_scoreCor = function(scoreList){
    count = 0
    corVal = NULL
    corEst = cor.test(as.numeric(scoreList[[1]]),as.numeric(scoreList[[2]]))
    corVal = c(corVal,c(corEst$estimate,corEst$p.value))
    return(corVal)
}

writeEMTscore<-function(gseID,GS76,KS){
    EMTDF=cbind(rownames(GS76),GS76,KS)
    EMTDF = as.data.frame(EMTDF)
    colnames(EMTDF) = c("Cellname","GS76_Score","KS_Score")
    fwrite(EMTDF,paste("../Output/",gseID,"_EMTScores.tsv",sep=''),sep='\t')
    }
