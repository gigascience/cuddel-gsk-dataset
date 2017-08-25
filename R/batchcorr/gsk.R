library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(specmine)
library(impute)
library(ggfortify)
library(batchCorr)

######################################
# Read in metadata for data analysis #
######################################

# Read metadata file
meta_all <- read.csv("/home/peter/gsk/meta/meta_all.csv")

# Use file_name_neg as rownames
rownames(meta_all) <- meta_all$file_name_neg
# Delete file_name_neg and file_name_pos columns
meta_all$file_name_neg <- NULL
meta_all$file_name_pos <- NULL

# Order files based on run order column
meta_all <- meta_all[order(meta_all$order),]

# Get paths for all negative files
neg_dir = "/home/peter/gsk/raw/esi_neg/netcdf"
neg_files <- character(0)
neg_files <- unlist(lapply(rownames(meta_all), function(x)
if (meta_all[x, "block"] == 1) {
    paste(neg_dir, "/block1neg/", x, ".cdf", sep = "")
}
else if (meta_all[x, "block"] == "2") {
    file <- paste(neg_dir, "/block2neg/", x, ".cdf", sep = "")
}
else if (meta_all[x, "block"] == "3") {
    file <- paste(neg_dir, "/block3neg/", x, ".cdf", sep = "")
}
else if (meta_all[x, "block"] == "4") {
    file <- paste(neg_dir, "/block4neg/", x, ".cdf", sep = "")
}
))

# Slice data
neg_qc_block1_filenames <- rownames(meta_all[which(meta_all$block == "1" & meta_all$type == "QC"),])
neg_qc_block2_filenames <- rownames(meta_all[which(meta_all$block == "2" & meta_all$type == "QC"),])
neg_qc_block3_filenames <- rownames(meta_all[which(meta_all$block == "3" & meta_all$type == "QC"),])
neg_qc_block4_filenames <- rownames(meta_all[which(meta_all$block == "3" & meta_all$type == "QC"),])

neg_sample_block1_filenames <- rownames(meta_all[which(meta_all$block == "1" & meta_all$type == "Sample"),])
neg_sample_block2_filenames <- rownames(meta_all[which(meta_all$block == "2" & meta_all$type == "Sample"),])
neg_sample_block3_filenames <- rownames(meta_all[which(meta_all$block == "3" & meta_all$type == "Sample"),])
neg_sample_block4_filenames <- rownames(meta_all[which(meta_all$block == "3" & meta_all$type == "Sample"),])

# Check
length(neg_files)
[1]371
all_files <- c(neg_qc_block1_filenames, neg_qc_block2_filenames, neg_qc_block3_filenames, neg_qc_block4_filenames)
all_files <- c(all_files, neg_sample_block1_filenames, neg_sample_block2_filenames, neg_sample_block3_filenames, neg_sample_block4_filenames)
length(all_files)
[1]371

# Create meta containing batch and grp information
sample_names <- substring(neg_files, 46)
sample_names <- gsub(".cdf", "", sample_names)
batch <- meta_all[sample_names, "block"]
batch <- paste("block", batch, sep = "")
grp <- meta_all[sample_names, "type"]
grp <- unlist(lapply(grp, as.character))
meta <- cbind(batch, grp)

########################
# Apply XCMS onto data #
########################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and will take time to complete!!
neg_xset <- xcmsSet(neg_files, step = 0.02, snthresh = 3, mzdiff = 0.05)

# Need to remove missing values in metadata otherwise XCMS group command fails
meta_all_no_na <- meta_all[, 1:11]
# Create batch and grp columns for use in batchcorr
meta_all_no_na$batch <- meta_all_no_na$block
meta_all_no_na$grp <- meta_all_no_na$type

# Add metadata to XCMS object
phenoData(neg_xset) <- meta_all_no_na

# Match peaks representing same analyte across samples
grp_neg_xset <- group(neg_xset, bw = 10, mzwid = 0.05)

# Perform retention time correction
# QC3=retcor(QC2, family="s", span=0.2)
retcor_neg_xset = retcor(grp_neg_xset, family = "s", span = 0.2)
retcor_neg_xsa <- xsAnnotate(retcor_neg_xset)
retcor_neg_peaklist <- getPeaklist(retcor_neg_xsa)

# After rt correction, the initial grouping is invalid so need to re-group
# QC_nofill=group(QC3,bw=1,mzwid=0.015,minfrac=.75)
nofill_neg_xset = group(retcor_neg_xset, bw = 10, mzwid = 0.05)
phenoData(nofill_neg_xset) <- meta
nofill_neg_xsa <- xsAnnotate(nofill_neg_xset)
nofill_neg_peaklist <- getPeaklist(nofill_neg_xsa)

# Fill in missing peaks
# QC_fill=fillPeaks(QC_nofill,method='chrom')
fill_neg_xset = fillPeaks(nofill_neg_xset, method = 'chrom')
fill_neg_xsa <- xsAnnotate(fill_neg_xset)
fill_neg_peaklist <- getPeaklist(fill_neg_xsa)

## Organise into peak table with missing data
# QCB=grabAlign(QC_nofill,batch='Batch_B',grp='QB')
# RefB=grabAlign(QC_nofill,batch='Batch_B',grp='Ref')
# QCF=grabAlign(QC_nofill,batch='Batch_F',grp='QF')
# RefF=grabAlign(QC_nofill,batch='Batch_F',grp='Ref')
# QCH=grabAlign(QC_nofill,batch='Batch_H',grp='QH')
# RefH=grabAlign(QC_nofill,batch='Batch_H',grp='Ref')
# PTnofill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
block1_QC <- grabAlign(nofill_neg_xset, batch='block1', grp='QC')
block1_sample <- grabAlign(nofill_neg_xset, batch='block1', grp='Sample')
block2_QC <- grabAlign(nofill_neg_xset, batch='block2', grp='QC')
block2_sample <- grabAlign(nofill_neg_xset, batch='block2', grp='Sample')
block3_QC <- grabAlign(nofill_neg_xset, batch='block3', grp='QC')
block3_sample <- grabAlign(nofill_neg_xset, batch='block3', grp='Sample')
block4_QC <- grabAlign(nofill_neg_xset, batch='block4', grp='QC')
block4_sample <- grabAlign(nofill_neg_xset, batch='block4', grp='Sample')
PTnofill <- rbind(block1_QC, block1_sample, block2_QC, block2_sample, block3_QC, block3_sample, block4_QC, block4_sample)

nofill_block1_qc = nofill_retcor_grp_neg_peaklist[, neg_qc_block1_filenames]
nofill_block1_samples = nofill_retcor_grp_neg_peaklist[, neg_sample_block1_filenames]
nofill_block2_qc = nofill_retcor_grp_neg_peaklist[, neg_qc_block2_filenames]
nofill_block2_samples = nofill_retcor_grp_neg_peaklist[, neg_sample_block2_filenames]
nofill_block3_qc = nofill_retcor_grp_neg_peaklist[, neg_qc_block3_filenames]
nofill_block3_samples = nofill_retcor_grp_neg_peaklist[, neg_sample_block3_filenames]
nofill_block4_qc = nofill_retcor_grp_neg_peaklist[, neg_qc_block4_filenames]
nofill_block4_samples = nofill_retcor_grp_neg_peaklist[, neg_sample_block4_filenames]
PTnofill = cbind(nofill_block1_qc, nofill_block1_samples, nofill_block2_qc, nofill_block2_samples, nofill_block3_qc, nofill_block3_samples, nofill_block4_qc, nofill_block4_samples)

## Organise into peaktable without missing data
# QCB=grabAlign(QC_fill,batch='Batch_B',grp='QB')
# RefB=grabAlign(QC_fill,batch='Batch_B',grp='Ref')
# QCF=grabAlign(QC_fill,batch='Batch_F',grp='QF')
# RefF=grabAlign(QC_fill,batch='Batch_F',grp='Ref')
# QCH=grabAlign(QC_fill,batch='Batch_H',grp='QH')
# RefH=grabAlign(QC_fill,batch='Batch_H',grp='Ref')
# PTfill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
fill_block1_qc = fill_retcor_grp_neg_peaklist[, neg_qc_block1_filenames]
fill_block1_samples = fill_retcor_grp_neg_peaklist[, neg_sample_block1_filenames]
fill_block2_qc = fill_retcor_grp_neg_peaklist[, neg_qc_block2_filenames]
fill_block2_samples = fill_retcor_grp_neg_peaklist[, neg_sample_block2_filenames]
fill_block3_qc = fill_retcor_grp_neg_peaklist[, neg_qc_block3_filenames]
fill_block3_samples = fill_retcor_grp_neg_peaklist[, neg_sample_block3_filenames]
fill_block4_qc = fill_retcor_grp_neg_peaklist[, neg_qc_block4_filenames]
fill_block4_samples = fill_retcor_grp_neg_peaklist[, neg_sample_block4_filenames]
PTfill = cbind(fill_block1_qc, fill_block1_samples, fill_block2_qc, fill_block2_samples, fill_block3_qc, fill_block3_samples, fill_block4_qc, fill_block4_samlsples)

## Perform batch alignment
# Extract peakinfo (i.e. m/z and rt of features)
peakIn <- cbind(nofill_retcor_grp_neg_peaklist[, "mz"], nofill_retcor_grp_neg_peaklist[, "rt"])
colnames(peakIn) <- c("mz", "rt")
rownames_peakIn <- integer(0)
for (i in 1:nrow(peakIn)) {
rownames_peakIn[i]<- paste("feature_", i, sep = "")
}
rownames(peakIn) <- rownames_peakIn

# Create meta object
batch <- c(rep('1', ncol(fill_block1_qc) + ncol(fill_block1_samples)), rep('2', ncol(fill_block2_qc) + ncol(fill_block2_samples)), rep('3', ncol(fill_block3_qc) + ncol(fill_block3_samples)), rep('4', ncol(fill_block4_qc) + ncol(fill_block4_samples)))
grp <- c(rep('Q', ncol(fill_block1_qc)), rep('S', ncol(fill_block1_samples)), rep('Q', ncol(fill_block2_qc)), rep('S', ncol(fill_block2_samples)), rep('Q', ncol(fill_block3_qc)), rep('S', ncol(fill_block3_samples)), rep('Q', ncol(fill_block4_qc)), rep('S', ncol(fill_block4_samples)))
meta <- cbind(batch, grp)

# Flag presence/missingness on batch level
bF <- batchFlag(t(PTnofill), meta, peakIn)

######################################################
# Find possible alignment candidates per sample type #
######################################################
# aIQ <- alignIndex2(bF, grpType = 'Q', mzdiff = 0.002, rtdiff = 15, report = T, reportName = 'splits_aIQ')
#alignIndex2=function(batchflag, flagType=c('Hard', 'Soft', 'All'), grpType='Q', mzdiff=0.005, rtdiff=10, report=TRUE, reportName='splits') {
grpType='Q'
mzdiff=0.002
rtdiff=15
report=T
reportName='splits_aIQ'

bF=batchflag
if (missing(flagType)) flagType='Hard'
if (flagType=='Hard') grpSub=bF$flagHard[bF$meta[, 2]==grpType,]# Take out matrix based on group type (such as Q for QC or R for Ref)
if (flagType=='Soft') grpSub=bF$flagSoft[bF$meta[, 2]==grpType,]# Take out matrix based on group type (such as Q for QC or R for Ref)
if (flagType=='All') grpSub=bF$flagAll[bF$meta[, 2]==grpType,]# Take out matrix based on group type (such as Q for QC or R for Ref)
mz=bF$mz
rt=bF$rt
head(grpSub)

a2=a1=align(grpSub, mz, rt, mzdiff=mzdiff,rtdiff=rtdiff)
a1Clust=a1$clusters
splits=which(a1Clust$dotProd!=0)
if (report==TRUE) pdf(file=paste(reportName, '.pdf', sep=''))
for (s in splits) {
cluster=a1$features[a1$features$cluster==s,]
# print(cluster)
mzSp=cluster$mz
rtSp=cluster$rt
flagSp=t(cluster[9:ncol(a1$features)])
alignSplit=clustSplit(flagSp, mzSp, rtSp)
text=paste('Original cluster', s)
if (report==TRUE) plotClust(bF, bF$meta[, 2]==grpType, cluster, text=text, color=alignSplit+1)
newClust=ifelse(alignSplit==1, s, max(a2$features$cluster)+alignSplit-1)
a2$features$cluster[a2$features$cluster==s]=newClust
}
if (report==TRUE) dev.off()
### Remove duplicates
a2$features=a2$features[a2$features$cluster%in%unique(a2$features$cluster[duplicated(a2$features$cluster)]),]
a2$features=a2$features[order(a2$features$cluster),]
a2$features$cluster=as.numeric(as.factor(a2$features$cluster))
a2$clusters=clusterMatrix(a2$features)
### Make index string for alignment
nClust=nrow(a2$clusters)
shiftList=numeric(ncol(bF$flagHard))
shiftGrp=character(ncol(bF$flagHard))
for (c in 1:nClust) {
feats=a2$features$featureIndex[a2$features$cluster==c]
shiftList[feats]=feats[1]
shiftGrp[feats]=as.character(grpType)
shift=data.frame(list=shiftList, grp=shiftGrp)
}
#return(list(grpType=grpType, shift=shift, events=a2$events, features=a2$features,clusters=a2$clusters, oldFeatures=a1$features, oldClusters=a1$clusters))
#}


##################################################


# Plot achieved alignments
plotAlign(bF, aIQ, plotType='pdf', reportName='clustPlots_aIQ')
# Perform alignment -> Peaktable
bA=batchAlign(bF, aIQ, PTfill, meta)
# Extract new peak table
PT=bA$PTalign

# Extract m/z and rt values




# Create CAMERA object
neg_xsa <- xsAnnotate(grp_neg_xset)
# Get peaklist
neg_peaklist <- getPeaklist(neg_xsa)
# Add numeric ID index into negative peaklist to keep track of peaks during data pre-treatment
neg_peaklist$idx <- seq.int(nrow(neg_peaklist))
# Move index to left hand side of data frame
neg_peaklist <- neg_peaklist[, c(ncol(neg_peaklist), 1:(ncol(neg_peaklist)-1))]
# Get list of neg QC filenames
get_neg_qc_file_paths <- function(x, output) {
file_name_neg <- x[3]
type <- x[10]
qc_meta_all <- x[12]
block <- x[13]
if (type == "QC" & block == "1") {
qc_file <- paste(neg_dir, "/block1neg/", file_name_neg, ".cdf", sep = "")
}
else if (type == "QC" & block == "2") {
qc_file <- paste(neg_dir, "/block2neg/", file_name_neg, ".cdf", sep = "")
}
else if (type == "QC" & block == "3") {
qc_file <- paste(neg_dir, "/block3neg/", file_name_neg, ".cdf", sep = "")
}
else if (type == "QC" & block == "4") {
qc_file <- paste(neg_dir, "/block4neg/", file_name_neg, ".cdf", sep = "")
}
else {
qc_file <- "Sample"
}
}
neg_qc_files <- apply(meta_all, 1, get_neg_qc_file_paths, output = qc_files)
# Remove "Sample" values in character vector
neg_qc_files <- neg_qc_files[! neg_qc_files %in% "Sample"]
# Sort negative QC filenames
neg_qc_files <- naturalsort(neg_qc_files)

# Combine neg_peaklist with sample_ID metadata from meta file
meta <- read.csv("/home/peter/gsk/meta/meta.csv")
neg_qc_fnames <- substr(neg_qc_files, 46, 80)
neg_qc_fnames <- substr(neg_qc_fnames, 1, nchar(neg_qc_fnames)-4)

########################
# Prepare data for PCA #
########################

# Get peak list df for neg QCs with index
# neg_qc_peaklist <- data.frame(neg_peaklist[1:12], neg_peaklist[neg_qc_fnames])
# Get peak list df for neg QCs
neg_qc_peaklist <- neg_peaklist[neg_qc_fnames]
###########################################
# Do PCA to check negative peak list data #
###########################################
# Create vector containing block information
sample_names <- colnames(neg_qc_peaklist)
block <- integer(0)
for (i in 1:length(sample_names)) {
if (grepl("block1", sample_names[i]) == 1) {
block[i]<- "block1"
}
else if (grepl("block2", sample_names[i]) == 1) {
block[i]<- "block2"
}
else if (grepl("block3", sample_names[i]) == 1) {
block[i]<- "block3"
}
else if (grepl("block4", sample_names[i]) == 1) {
block[i]<- "block4"
}
}
# Transpose data
pca_data <- t(neg_qc_peaklist)
# Add block information to PCA data
pca_data <- cbind(pca_data, block)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[, colSums(is.na(pca_data)) == 0]
write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[, 1:ncol(pca_data)-1]), data = pca_data, colour = 'block', main = 'PCA on unprocessed negative QC data')
ggsave("unprocessed_neg_qc_pca.png")






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
data(PTfill, PTnofill, meta, BData, FData,HData)
##########################
## Perform batch alignment
# Extract peakinfo (i.e. m/z and rt of features)
peakIn=peakInfo(PTnofill)
# Flag presence/missingness on batch level
bF=batchFlag(PTnofill, meta, peakIn)
# Find possible alignment candidates per sample type
aIQ=alignIndex(bF, grpType='Q', mzdiff=0.002,rtdiff=15, report=T, reportName='splits_aIQ')
# Plot achieved alignments
plotAlign(bF, aIQ, plotType='pdf', reportName='clustPlots_aIQ')
# Perform alignment -> Peaktable
bA=batchAlign(bF, aIQ, PTfill, meta)
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
B=driftWrap(BData, refType='one', report=TRUE)
F=driftWrap(FData, refType='one', report=TRUE)
H=driftWrap(HData, refType='one', report=TRUE)

library(VennDiagram)
venn.diagram(list(B=B$finalVars, F=F$finalVars, H=H$finalVars), file='VennBatch.png')

png(file='histVars.png', width=1024, height=1024, pointsize = 24)
par(mar=c(4, 4, 0, 0)+0.2)
hist(cv(B$QCFeats), 200, col=rgb(0.1, 0.1, 0.1, 0.5), xlim=c(0, 1), main='', xlab='CV (feature)')
hist(cv(B$QCFeatsCorr), 200, col=rgb(0.9, 0.9, 0.9, .6), add=TRUE)
# hist(cv(B$QCFeatsFinal),20,col=rgb(0.3,0.3,0.3,.3),add=TRUE)
lines(rep(0.3, 2), c(0, 800), lty=2)
text(0.3, 800, 'CV limit 30%', pos=4)
box(bty='l')
legend('topright', legend=c('Before correction', 'After correction'), fill=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.9, 0.9, 0.9, .6)))
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

batchObjs=list(B, F, H)
PTComb=batchComb(batchObjs)  # Combine features present in all batches

#####     ---------::::::::::::: Run to here


################################################
## Perform between batch intensity normalisation
## By long-term reference samples

refs=refOut(PTComb, meta) # Aggregate ref sample info on batch level
RC=refCorr(PTComb, meta, refs)  # Perform between batch normalisation for those features passing Ref heuristic criterion
PT_RefNorm=RC$PTRef  # Extract reference normalised peak table

# NB: Final normalization is not complete!!!
# Batches where Ref heuristic criterion is not passed still need to be population normalized...
# Info on which batches need to be population normalized is found in RC$refCorr:
# TRUE denotes normalisation by reference samples
# FALSE denotes NOT normalised -> Need to be population normalised instead.
#
# This is not performed for this data set, since sample populations are not present,
# but only batch-specific QC and long-term reference samples.
