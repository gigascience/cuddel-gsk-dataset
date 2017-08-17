###################
## Read in raw data
#   The indented section was already performed to generate data available in package
    # filepath="C:/R/QCData/OrgData"  # Specify directory of centroid files
    ## Read samples
    # QC1=xcmsSet(filepath,method='centWave',prefilter=c(3,440),peakwidth=c(5,76),snthresh=6,mzdiff=0.0045,ppm=15)
    # QC2=group(QC1,bw=10,minfrac=0.75,minsamp=1,mzwid=0.015,sleep=0)
    # QC3=retcor(QC2, family="s", span=0.2)
    # QC_nofill=group(QC3,bw=1,mzwid=0.015,minfrac=.75)
    # QC_fill=fillPeaks(QC_nofill,method='chrom')
    ## Organise into peaktable with missing data
    # QCB=grabAlign(QC_nofill,batch='Batch_B',grp='QB')
    # RefB=grabAlign(QC_nofill,batch='Batch_B',grp='Ref')
    # QCF=grabAlign(QC_nofill,batch='Batch_F',grp='QF')
    # RefF=grabAlign(QC_nofill,batch='Batch_F',grp='Ref')
    # QCH=grabAlign(QC_nofill,batch='Batch_H',grp='QH')
    # RefH=grabAlign(QC_nofill,batch='Batch_H',grp='Ref')
    # PTnofill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
    ## Organise into peaktable without missing data
    # QCB=grabAlign(QC_fill,batch='Batch_B',grp='QB')
    # RefB=grabAlign(QC_fill,batch='Batch_B',grp='Ref')
    # QCF=grabAlign(QC_fill,batch='Batch_F',grp='QF')
    # RefF=grabAlign(QC_fill,batch='Batch_F',grp='Ref')
    # QCH=grabAlign(QC_fill,batch='Batch_H',grp='QH')
    # RefH=grabAlign(QC_fill,batch='Batch_H',grp='Ref')
    # PTfill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
    ## Set up metadata (Quick'n'Dirty approach)
    # batch=c(rep('B',nrow(QCB)+nrow(RefB)),rep('F',nrow(QCF)+nrow(RefF)),rep('H',nrow(QCH)+nrow(RefH)))
    # grp=c(rep('Q',nrow(QCB)),rep('R',nrow(RefB)),rep('Q',nrow(QCF)),rep('R',nrow(RefF)),rep('Q',nrow(QCH)),rep('R',nrow(RefH)))
    # meta=cbind(batch,grp)
    ## Save relevant data files
    # save(PTfill,file='PT_fill.RData')
    # save(PTnofill,file='PT_nofill.RData')
    # save(meta,file='meta.RData')

# LOAD DATA
library(batchCorr)
data(PTfill,PTnofill,meta,BData,FData,HData)
##########################
## Perform batch alignment
# Extract peakinfo (i.e. m/z and rt of features)
peakIn=peakInfo(PTnofill)
# Flag presence/missingness on batch level
bF=batchFlag(PTnofill,meta,peakIn)
# Find possible alignment candidates per sample type
aIQ=alignIndex(bF,grpType='Q',mzdiff=0.002,rtdiff=15,report=T,reportName='splits_aIQ')
# Plot achieved alignments
plotAlign(bF,aIQ,plotType='pdf',reportName='clustPlots_aIQ')
# Perform alignment -> Peaktable
bA=batchAlign(bF,aIQ,PTfill,meta)
# Extract new peak table
PT=bA$PTalign

dim(PT)
# [1]    90 11284









#########################################
## Perform within batch drift corrections
#   The indented section was already performed to generate data available in package
    # Extract data per batch
    # BData=grabWrapBA(QC_fill,PT,batch='Batch_B',QC='QB')
    # FData=grabWrapBA(QC_fill,PT,batch='Batch_F',QC='QF')
    # HData=grabWrapBA(QC_fill,PT,batch='Batch_H',QC='QH')
# Perform drift correction within batch
B=driftWrap(BData,refType='one',report=TRUE)
F=driftWrap(FData,refType='one',report=TRUE)
H=driftWrap(HData,refType='one',report=TRUE)

library(VennDiagram)
venn.diagram(list(B=B$finalVars,F=F$finalVars,H=H$finalVars),file='VennBatch.png')

png(file='histVars.png',width=1024,height=1024,pointsize = 24)
par(mar=c(4,4,0,0)+0.2)
hist(cv(B$QCFeats),200,col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,1),main='',xlab='CV (feature)')
hist(cv(B$QCFeatsCorr),200,col=rgb(0.9,0.9,0.9,.6),add=TRUE)
# hist(cv(B$QCFeatsFinal),20,col=rgb(0.3,0.3,0.3,.3),add=TRUE)
lines(rep(0.3,2),c(0,800),lty=2)
text(0.3,800,'CV limit 30%',pos=4)
box(bty='l')
legend('topright',legend=c('Before correction','After correction'),fill=c(rgb(0.1,0.1,0.1,0.5),rgb(0.9,0.9,0.9,.6)))
dev.off()

## save(B,F,H,file='batchDrift.Rdata')
sum(cv(B$QCFeats)<0.3)
sum(cv(B$QCFeatsCorr)<0.3)
mean(cv(B$QCFeats))
mean(cv(B$QCFeatsCorr))
mean(cv(B$QCFeatsFinal))

sum(cv(F$QCFeats)<0.3)
sum(cv(F$QCFeatsCorr)<0.3)
mean(cv(F$QCFeats))
mean(cv(F$QCFeatsCorr))
mean(cv(F$QCFeatsFinal))

sum(cv(H$QCFeats)<0.3)
sum(cv(H$QCFeatsCorr)<0.3)
mean(cv(H$QCFeats))
mean(cv(H$QCFeatsCorr))
mean(cv(H$QCFeatsFinal))

##########################################################
## Combine batch data and restrict to common features only

batchObjs=list(B,F,H)
PTComb=batchComb(batchObjs)  # Combine features present in all batches

#####     ---------::::::::::::: Run to here


################################################
## Perform between batch intensity normalisation
## By long-term reference samples

refs=refOut(PTComb,meta) # Aggregate ref sample info on batch level
RC=refCorr(PTComb,meta,refs)  # Perform between batch normalisation for those features passing Ref heuristic criterion
PT_RefNorm=RC$PTRef  # Extract reference normalised peak table

# NB: Final normalization is not complete!!!
# Batches where Ref heuristic criterion is not passed still need to be population normalized...
# Info on which batches need to be population normalized is found in RC$refCorr:
# TRUE denotes normalisation by reference samples
# FALSE denotes NOT normalised -> Need to be population normalised instead.
#
# This is not performed for this data set, since sample populations are not present,
# but only batch-specific QC and long-term reference samples.
