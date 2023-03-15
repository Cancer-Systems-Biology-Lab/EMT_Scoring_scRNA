##==================================================================================================
## Raw read count to TPM conversion code
## Author: Priyanka Chakraborty
## Date: 06-08-2020
## References : 
##  https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
##  http://robpatro.com/blog/?p=235
##  https://www.biostars.org/p/390038/
##====================================================================================================



countToTpm <- function(rawCount, GSEID) {
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
    # Get the common gene indices for gene length annotation
    idx2 = match(commonGenes,geneLengthMap[,2])
    # Subset gene length annotation for common genes
    geneLengthSubset = geneLengthMap[idx2, ]
    geneSymbol = geneLengthSubset[ ,geneCol]
    featureLength = geneLengthSubset[ ,3]

    rownames(rawCount) = geneSymbol
    Samples = colnames(rawCount) 
    # Ensure valid arguments.
    stopifnot(length(featureLength) == nrow(rawCount))

    # Compute effective lengths of features in each library.
    effLen <- featureLength

    # Process one column at a time.
    tpm <- do.call(cbind, lapply(1:ncol(rawCount), function(i) {
        rate = (rawCount[,i])/effLen
        rate/sum(rate, na.rm = T) * 1e6
    }))

    tpm = log2(tpm+1)
    # Copy the row and column names from the original matrix.
    colnames(tpm) <- Samples 
    rownames(tpm) <- toupper(geneSymbol)

    outFile = paste("../Data_generated/", GSEID, "_TPM.tsv", sep = "")
    write.table(tpm, outFile, sep = '\t', quote = F)

    return(tpm)

}

rnaToMA = function(TPMVal){
    MA = (0.57 +(TPMVal* 0.37 ))
    return(MA)
}