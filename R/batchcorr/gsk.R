library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(specmine)
library(impute)
library(ggfortify)
library(batchCorr)

# To see how to create a phenodata object
library(Biobase)
data(sample.ExpressionSet)
head(pData(sample.ExpressionSet))
     sex    type score  newColumn
A Female Control  0.75 -0.3141766
B   Male    Case  0.40  2.3844318
C   Male Control  0.73  1.0936301
D   Male    Case  0.42 -2.6438696
E Female    Case  0.93 -0.3369824
F   Male Control  0.22  0.6235963

pData(sample.ExpressionSet)$newColumn = rnorm(ncol(sample.ExpressionSet))
head(pData(sample.ExpressionSet))
> colnames(exprs(sample.ExpressionSet))
 [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S"
[20] "T" "U" "V" "W" "X" "Y" "Z"
> head(exprs(sample.ExpressionSet))
> head(exprs(sample.ExpressionSet))
                       A         B        C        D        E       F        G
AFFX-MurIL2_at  192.7420  85.75330 176.7570 135.5750 64.49390 76.3569 160.5050
AFFX-MurIL10_at  97.1370 126.19600  77.9216  93.3713 24.39860 85.5088  98.9086
AFFX-MurIL4_at   45.8192   8.83135  33.0632  28.7072  5.94492 28.2925  30.9694
AFFX-MurFAS_at   22.5445   3.60093  14.6883  12.3397 36.86630 11.2568  23.0034
AFFX-BioB-5_at   96.7875  30.43800  46.1271  70.9319 56.17440 42.6756  86.5156
AFFX-BioB-M_at   89.0730  25.84610  57.2033  69.9766 49.58220 26.1262  75.0083

# Therefore need to create a phenodata object like this
					batch	grp
GSK_neg_block1_09r
GSK_neg_block1_10r
GSK_neg_block1_16r
GSK_neg_block1_21r
GSK_neg_block1_26r
GSK_neg_block1_31r


# Then need to add DF when creating xcmsSet object i.e. xcmsSet(neg_files, step = 0.02, snthresh=3, mzdiff = 0.05)

# Need to consider batch differences which might affect peaks in QCs due to the
# 4 analytical blocks are removed by XCMS during peak alignment.

######################################
# Read in metadata for data analysis #
######################################

# Read metadata file
meta_all <- read.csv("/home/peter/gsk/meta/meta_all.csv")

# Use file_name_neg as rownames
rownames(meta_all) <- meta_all$file_name_neg

# Delete file_name_neg and file_name_pos columns
meta_all$file_name_neg <- NULL
meta_all$file_name_pos <-

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

########################
# Apply XCMS onto data #
########################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and will take time to complete!!
neg_xset <- xcmsSet(neg_files, step = 0.02, snthresh = 3, mzdiff = 0.05)

# Need to remove missing values in metadata otherwise XCMS group cmd fails
meta_all_no_na <- meta_all[, 1:11]
# Add metadata to XCMS object
phenoData(neg_xset) <- meta_all_no_na

# Match peaks representing same analyte across samples
grp_neg_xset <- group(neg_xset, bw = 10, mzwid = 0.05)

# QC3=retcor(QC2, family="s", span=0.2)
retcor_grp_neg_xset = retcor(grp_neg_xset, family = "s", span = 0.2)

# QC_nofill=group(QC3,bw=1,mzwid=0.015,minfrac=.75)
# nofill_retcor_grp_neg_xset = group(retcor_grp_neg_xset, bw = 1, mzwid = 0.015, minfrac = .75)
nofill_retcor_grp_neg_xset = group(retcor_grp_neg_xset, bw = 10, mzwid = 0.05)

# QC_fill=fillPeaks(QC_nofill,method='chrom')
fill_retcor_grp_neg_xset = fillPeaks(nofill_retcor_grp_neg_xset, method = 'chrom')

## Organise into peak table with missing data
# QCB=grabAlign(QC_nofill,batch='Batch_B',grp='QB')
neg_qc_block1 = grabAlign(nofill_retcor_grp_neg_xset, batch = 'QC_neg_block1,', grp = 'block1neg')

# QCF=grabAlign(QC_nofill,batch='Batch_F',grp='QF')
# QCH=grabAlign(QC_nofill,batch='Batch_H',grp='QH')
# PTnofill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)

# Create CAMERA object
neg_xsa <- xsAnnotate(grp_neg_xset)
# Get peaklist
neg_peaklist <- getPeaklist(neg_xsa)
# Add numeric ID index into negative peaklist to keep track of peaks during data pre-treatment
neg_peaklist$idx <- seq.int(nrow(neg_peaklist))
# Move index to left hand side of data frame
neg_peaklist <- neg_peaklist[,c(ncol(neg_peaklist),1:(ncol(neg_peaklist)-1))]

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
neg_qc_files <- neg_qc_files [! neg_qc_files %in% "Sample"]
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
    block[i] <- "block1"
  }
  else if (grepl("block2", sample_names[i]) == 1) {
    block[i] <- "block2"
  }
  else if (grepl("block3", sample_names[i]) == 1) {
    block[i] <- "block3"
  }
  else if (grepl("block4", sample_names[i]) == 1) {
    block[i] <- "block4"
  }
}

# Transpose data
pca_data <- t(neg_qc_peaklist)
# Add block information to PCA data
pca_data <- cbind(pca_data, block)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = 'block', main = 'PCA on unprocessed negative QC data')
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
