# EMT_Scoring_RNASeq

This code takes raw read counts from the user and calculated TPM expression data and calculates EMT score using methods mentioned in Byers et al.,2013, Tan et al., 2014, George et al., 2017 and Chakraborty., et al 2020

### Preprocessing of the raw data

After QC filtering, align using STAR-aligner using appropriate reference genome (hg38/mm10) and then calculate raw read counts using htseq-count. 
Alternatively, you can download the raw counts from GEO.

### Requirements
- R:
		- dplyr
		- data.table
		- readxl
		- matlabr
- MATLAB (If you do not have MATLAB and MLRScore is not required, checkout the NoMLR branch which does not have this dependency)

### Running the code

Please run the following in a terminal :

```
Rscript  ./all_GSE_EMT.R
```

This code will calculate the TPM expression data, EMT scores of three metrics: 76GS, KS and MLR. Also the correlations between these three metrics.

***If you have any problem running this code, please let me know***

