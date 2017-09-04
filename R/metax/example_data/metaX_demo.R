options(bitmapType="cairo");
#load metaX package
library(metaX)
para <- new("metaXpara")
#input the path of peak intensity file
pfile <- "./Input_data.csv"
#input the path of sample list file which containing the sample,group,batch and order information
sfile <- "./sample_list.txt"
#set the outdir name
para@outdir <- "metaX_result_neg"
#mode
para@prefix <- "neg"
para@sampleListFile <- sfile
#comparation pairs
para@ratioPairs <-"A:B;"
para <- importDataFromXCMS(para,file=pfile)
plsdaPara <- new("plsDAPara")
#set the scaling method for plsda
plsdaPara@scale = "pareto"
plsdaPara@cpu = 4
res <- doQCRLSC(para,cpu=1)
p <- metaXpipe(para,plsdaPara=plsdaPara,missValueRatioQC = 0.5, missValueRatioSample = 0.8,cvFilter=0.3,idres=NULL,qcsc=1,scale="pareto",remveOutlier=FALSE,nor.method="pqn",t=1,nor.order = 1,pclean = FALSE,doROC=FALSE) #set the other parameters
save(p,file="p1.rda")
sessionInfo()

