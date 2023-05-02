library(trqwe)
library(Seurat)
library(data.table)

countToTpm_sc <- function(Sample) {
    # Read count matrix
    rawCount<-mcreadRDS(paste0(Sample,"_counts.rds"), mc.cores=4)
    rawCount <- cbind(rownames(rawCount), rawCount)
    # Read gene length annotation
    geneLengthMap_hg = read.table(file="../Annotation/hg38_gene_length_kb.txt", sep = "\t",header = F,stringsAsFactors = F)
    commonGenes_hg = intersect(rawCount[,1], geneLengthMap_hg[,2])
    geneLengthMap_mm = read.table(file="../Annotation/mm10_gene_length_kb.txt", sep = "\t",header = F,stringsAsFactors = F)
    commonGenes_mm = intersect(rawCount[,1], geneLengthMap_mm[,2])
    geneCol = 2

    if(length(commonGenes_hg) > length(commonGenes_mm)) {
        geneLengthMap = geneLengthMap_hg
        commonGenes = commonGenes_hg
    }
    else{
        geneLengthMap = geneLengthMap_mm
        commonGenes = commonGenes_mm
    }
    # Get the common gene indices for count matrix
    idx1 = match(commonGenes,rawCount[,1])
    # Subset count matrix for common genes
    rawCount = rawCount[idx1,-c(1)]
    gc()
    # Get the common gene indices for gene length annotation
    idx2 = match(commonGenes,geneLengthMap[,2])
    # Subset gene length annotation for common genes
    geneLengthSubset = geneLengthMap[idx2, ]
    geneSymbol = geneLengthSubset[ ,geneCol]
    featureLength = geneLengthSubset[ ,3]

    rownames(rawCount) = geneSymbol
    CellNames = colnames(rawCount) 
    # Ensure valid arguments.
    stopifnot(length(featureLength) == nrow(rawCount))

    # Compute effective lengths of features in each library.
    effLen <- featureLength

    # Process one column at a time.
    tpm <- do.call(cbind, lapply(1:ncol(rawCount), function(i) {
    rate = (rawCount[,i])/effLen
    rate/sum(rate, na.rm = T) * 1e6
    }))

    rm(rawCount)
    gc()

    tpm = log2(tpm+1)
    gc()
    # Copy the row and column names from the original matrix.
    colnames(tpm) <- CellNames 
    rownames(tpm) <- toupper(geneSymbol)
    # Set na values to 0
    tpm[is.na(tpm)] = 0
    # Save TPM
    mcsaveRDS(tpm, file = paste0('../Data_generated/',Sample,'_TPM.rds'), mc.cores = 4)
}

rnaToMA = function(TPMVal){
    MA = (0.57 +(TPMVal* 0.37 ))
    return(MA)
}