# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 4/9/2017

library("metaX")

# Load example peak intensity data from MTBLS79
para <- new("metaXpara")
pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
head(para@rawPeaks[,1:4])

# Set the sample list file path
sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
sampleListFile(para) <- sfile

# Run integrated analysis
ratioPairs(para) <- "S:C"
plsdaPara <- new("plsDAPara")
para@outdir <- "output"
p <- metaXpipe(para,plsdaPara=plsdaPara,cvFilter=0.3,remveOutlier = TRUE, outTol = 1.2, doQA = TRUE, doROC = TRUE, qcsc = 1, nor.method = "pqn", pclean = TRUE, t = 1, scale = "uv",)

# Normalisation to correct data within batch experiment analytical variation, and batch-to-batch variation
> para <- new("metaXpara")
> pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
> file <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
> rawPeaks(para) <- read.delim(pfile,check.names = FALSE)[1:100,]
> sampleListFile(para) <- sfile
> para <- reSetPeaksData(para)
> para <- missingValueImpute(para)
missingValueImpute: value
Mon Sep  4 15:03:56 2017        Missing value imputation for 'value'
Missing value in total: 781
Missing value in QC sample: 139
Missing value in non-QC sample: 642
Mon Sep  4 15:03:56 2017        The ratio of missing value: 4.5407%
<=0: 0
Missing value in total after missing value inputation: 0
<=0 value in total after missing value inputation: 0

> res <- doQCRLSC(para,cpu=1)

