
#################################
# Use metaXpipe to analyse data #
#################################

para <- new("metaXpara")
pfile <- "/home/peter/gsk/raw/esi_neg/netcdf/test.tab"
sfile <- "/home/peter/gsk/raw/esi_neg/netcdf/sampleListFile.tab"
rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
sampleListFile(para) <- sfile
ratioPairs(para) <- "A:B;C:D;A:C;B:D"
plsdaPara <- new("plsDAPara")
para@outdir <- "output"
p <- metaXpipe(para,plsdaPara=plsdaPara,cvFilter=0.3,remveOutlier = TRUE, outTol = 1.2, doQA = TRUE, doROC = TRUE, qcsc = 1, nor.method = "pqn", pclean = TRUE, t = 1, scale = "uv",)


