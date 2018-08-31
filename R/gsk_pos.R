# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 30/7/2018

library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(mice)
library(ggfortify)
library(metaX)
library(VIM)
library(ropls)
library(MAIT)
library(Morpho)

# source("functions.R")

# Location of GSK data set
datadir = "/home/peter/"
# File path for positive files
pos_dir = paste(datadir, "gsk/raw/esi_pos/netcdf", sep="")

# Output path
output_path <- paste(pos_dir, "/output", sep="")

######################################
# Read in metadata for data analysis #
######################################

# Read in and sort metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta_20180831.csv", sep=""))
meta <- meta[naturalorder(meta$file_name_pos),]
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all_20180831.csv", sep=""))
meta_all <- meta_all[naturalorder(meta_all$file_name_pos), ]

# Get list of positive file paths
pos_file_paths <- getAllGSKFilePaths(mode="positive")

# Get list of positive files
pos_files <- meta[, "file_name_pos"]
pos_files <- as.character(pos_files)

# There are 371 positive files
length(pos_files)
## [1] 371

# 21 of the 371 positive samples are duplicated
dup <- pos_files[duplicated(pos_files)]
# All these 21 samples are QC samples
## > getQCSampleMetadata(dup, mode="positive")
## [1] "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC" "QC"

# Remove duplicated file paths for XCMS processing
pos_file_paths <- unique(pos_file_paths)

# Get unique positive sample names
unique_pos_files <- unique(pos_files)


######################################
# Apply XCMS onto QC and sample data #
######################################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and takes time to complete!!
pos_xset <- xcmsSet(pos_file_paths, step=0.02, snthresh=3, mzdiff=0.05)

# Match peaks representing same analyte across samples
grp_pos_xset <- group(pos_xset, bw=10, mzwid=0.05)

# Create CAMERA object
pos_xsa <- xsAnnotate(grp_pos_xset)
# Get peaklist
pos_peaklist <- getPeaklist(pos_xsa)
# Add numeric index to positive peaklist to keep track of peaks during data pre-treatment
pos_peaklist$idx <- seq.int(nrow(pos_peaklist))
# Move index to left hand side of data frame
pos_peaklist <- pos_peaklist[, c(ncol(pos_peaklist), 1:(ncol(pos_peaklist)-1))]

# Output pos_peaklist data
write.table(
    pos_peaklist,
    file = paste(output_path, "/xcms_pos_peaklist.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)

#############################
# Load saved peak list data #
#############################

xcms_pos_peaks = read.csv(paste(output_path, "/xcms_pos_peaklist.csv", sep=""))


#####################################################
# Do PCA on unprocessed positive block 3 QC samples #
# which have been used as block 2 QC samples and    #
# block 2 plasma samples                            #
#####################################################

## length(dup)
## [1] 21

# Get block 2 plasma sample names
block2_plasma_sample_names <- getPlasmaSampleNamesByBlock(mode='positive', block='2')


# Prepare pos_peaklist data
pca_data <- cbind(xcms_pos_peaks[, dup], xcms_pos_peaks[, block2_plasma_sample_names])
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
sample_meta <- getQCSampleMetadata(samples, mode="positive")
pca_meta <- cbind(block_meta, sample_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_meta',
    main='PCA on unprocessed positive block2 sample data and its block3 QC samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_pos_qc_sample_block2_samples_pca.png", sep=""))

autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on unprocessed positive block2 sample data and its block3 QC samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_pos_qc_sample_block2_samples_pca.png", sep=""))


#################################################
# Do PCA on unprocessed positive QC and samples #
#################################################

# Prepare pos_peaklist data
pca_data <- xcms_pos_peaks[, unique_pos_files]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
sample_meta <- getQCSampleMetadata(samples, mode="positive")
pca_meta <- cbind(block_meta, sample_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_meta',
    main='PCA on unprocessed positive QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_pos_qc_sample_pca_shape_block_colour_sample.png", sep=""))

autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on unprocessed positive QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_pos_qc_sample_pca.png", sep=""))

# There is an outlier which is affect the PCA display
# Need to remove this outlier


##############################
# Do PCA on positive QC data #
##############################

# Prepare QC positive peaklist data
meta_all_qc_rows <- meta_all[, "type"]=="QC"
pos_qc_names <- meta_all[meta_all_qc_rows, "file_name_pos"]
pos_qc_names <- as.character(pos_qc_names)
# Subset QC data from positive peak list
pca_data <- xcms_pos_peaks[, pos_qc_names]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on unprocessed positive QC samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_pos_qc_pca.png", sep=""))


############################################################################
# Calculate percentage of missing values for each feature row for QCs only #
############################################################################

# Prepare QC positive data
qc_pos_peaklist <- xcms_pos_peaks[, pos_qc_names]
# Calculate percentage of missing values by peak row
percent_nas <- rowMeans(is.na(qc_pos_peaklist))
# Add columns
idx <- rownames(qc_pos_peaklist)
qc_pos_peaklist <- cbind(idx, percent_nas, qc_pos_peaklist)
# Remove rows containing 40% or more missing values
na_filtered_qc_pos_peaklist <- qc_pos_peaklist[qc_pos_peaklist$percent_nas < .40, ]

# qc_pos_peaklist has 6198 rows. After filtering rows with 40% na, na_filtered_qc_pos_peaklist contains 4574 rows


###############################################################
# Do PCA after removing features with 40% missing data values #
###############################################################

# Prepare positive peaklist data
nr <- colnames(na_filtered_qc_pos_peaklist) %in% c("percent_nas", "idx")
pca_data <- na_filtered_qc_pos_peaklist[, !nr]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on positive QCs after peaks missing 40% values removed',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/40_percent_filtered_pos_qc_pca.png", sep=""))


###################################################################
# Impute k-nearest neighbours to replace remaining missing values #
###################################################################

# Prepare data for knn analysis
nr <- colnames(na_filtered_qc_pos_peaklist) %in% c("percent_nas", "idx")
impute_data <- na_filtered_qc_pos_peaklist[, !nr]
# Transpose data
impute_data <- t(impute_data)

# Identify which column peak features have missing values
na_peak_cols <- colnames(impute_data)[colSums(is.na(impute_data)) > 0]
# Perform K-means clustering using VIM - takes time!!
knn <- kNN(impute_data, variable=na_peak_cols)
# Copy peak features from kNN result
knn_qc_pos_peaklist <- knn[, 1:ncol(impute_data)]
knn_qc_pos_peaklist <- as.matrix(knn_qc_pos_peaklist)
# Copy rownames and column names from impute_knn_data to rownames for knn
rownames(knn_qc_pos_peaklist) <- rownames(impute_data)
colnames(knn_qc_pos_peaklist) <- colnames(impute_data)
# Transpose
knn_qc_pos_peaklist <- t(knn_qc_pos_peaklist)


#####################################################
# Do PCA plot to check effect of imputed knn values #
#####################################################
pca_data <- knn_qc_pos_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on positive QC data with imputed k-means for missing values',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/impute_knn_pos_qc_pca.png", sep=""))


######################################
# Filter k-means imputed data by RSD #
######################################

# Prepare data for RSD analysis
rsd_qc_pos_peaklist <- knn_qc_pos_peaklist

# Calculate RSD
rowRSD <- apply(rsd_qc_pos_peaklist, 1, function(data) {
    sd <- sd(data)
    mean <- mean(data)
    rsd <- (sd/mean) * 100
    return(rsd)
})

# Find out which rows are over 20% RSD and remove them
bad_rows <- which(rowRSD > 20)
rsd_qc_pos_peaklist <- rsd_qc_pos_peaklist[-bad_rows, ]


###########################################
# Do PCA to check effect of RSD filtering #
###########################################

pca_data <- rsd_qc_pos_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on RSD-filtered QC positive data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/rsd_pos_qc_pca.png", sep=""))


##########################################################################
# Combine feature vectors from plasma samples with filtered positive QCs #
##########################################################################

# Get peak indices of pretreated QC feature vectors
pretreated_qc_idx <- rownames(rsd_qc_pos_peaklist)
# Get plasma peaklist from pos_peaklist based on pretreated QC indices
filtered_pos_peaklist <- pos_peaklist[rownames(pos_peaklist) %in% pretreated_qc_idx, ]
pos_plasma_sample_names <- getPlasmaSampleNames(mode="positive")
filtered_plasma_peaklist <- filtered_pos_peaklist[, colnames(filtered_pos_peaklist) %in% pos_plasma_sample_names]
# Combine QC and plasma samples from positive peaklist
filtered_qc_plasma_peaklist <- cbind(rsd_qc_pos_peaklist, filtered_plasma_peaklist)


####################################################
# Do PCA on filtered QC and plasma sample peaklist #
####################################################

pca_data <- filtered_qc_plasma_peaklist
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Sample type information
sample_type_meta <- getQCSampleMetadata(samples, mode="positive")
# Convert block_meta character vector to matrix
pca_meta <- cbind(block_meta, sample_type_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='sample_type_meta',
    colour='block_meta',
    main='PCA on filtered positive QCs with positive plasma samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/filtered_pos_qc_plasma_samples_pca.png", sep=""))

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_type_meta',
    main='PCA on pretreated positive QCs with positive plasma samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/filtered_pos_qc_plasma_samples_pca_shape_block_colour_sample.png", sep=""))


###########################################################
# Use k-means to impute missing values for plasma samples #
###########################################################

# Prepare data for imputation
impute_data <- filtered_qc_plasma_peaklist
impute_data <- t(impute_data)

# Identify column peak features with missing values
na_cols <- colnames(impute_data)[colSums(is.na(impute_data)) > 0]
# Perform K-means clustering using VIM
knn <- kNN(impute_data, variable=na_cols)

# Copy peak features from kNN result
knn_qc_plasma_peaklist <- knn[, 1:ncol(impute_data)]
knn_qc_plasma_peaklist <- as.matrix(knn_qc_plasma_peaklist)
# Copy rownames and column names
rownames(knn_qc_plasma_peaklist) <- rownames(impute_data)
colnames(knn_qc_plasma_peaklist) <- colnames(impute_data)

# Transpose
knn_qc_plasma_peaklist <- t(knn_qc_plasma_peaklist)


#########################################################################
# Do PCA on combined rsd_qc_neg_peaklist and pretreated_sample_peaklist #
#########################################################################

pca_data <- knn_qc_plasma_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Sample type information
sample_type_meta <- getQCSampleMetadata(samples, mode="positive")
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta, sample_type_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_type_meta',
    main='PCA on pretreated positive QCs with knn positive samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pretreated_qc_knn_pos_samples_pca1.png", sep=""))

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='sample_type_meta',
    colour='block_meta',
    main='PCA on pretreated positive QCs with knn positive samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pretreated_qc_knn_pos_samples_pca2.png", sep=""))


############################################################
# Use metaX data normalisation to correct within batch and #
# batch-to-batch variation in positive QCs and samples     #
############################################################

# Need to consider batch differences affecting peaks in QCs
# due to the 4 analytical blocks are removed by XCMS during
# peak alignment.

# Prepare data
signal_corr_data <- knn_qc_plasma_peaklist
signal_corr_data <- cbind(rownames(signal_corr_data), signal_corr_data)
colnames(signal_corr_data)[1] <- "name"
signal_corr_data <- data.frame(signal_corr_data)

# Some sample names contain ".1" which might be causing problems
# for the normalisation step. These ".1" sample names are in
# signal_corr_data and can be removed as follows:

signal_corr_data_colnames <- colnames(signal_corr_data)
bad_cols <- which(grepl("\\.1", signal_corr_data_colnames))
signal_corr_data <- signal_corr_data[,-bad_cols]

# Check
head(signal_corr_data[, 1:4])
##     name GSK_pos_block1_81 GSK_pos_block1_86 GSK_pos_block3_82
## 20    20  77943.6834000002  77943.6834000002  77943.6834000002
## 37    37  89185.3310418605  89185.3310418605  89185.3310418605
## 56    56  1409220.34612157  1409220.34612157  1485652.15393913
## 68    68  321244.606438095  756865.658822728  761105.699639999
## 70    70  51248.4487704545  51248.4487704545  60061.8399727272
## 105  105  249396.199392857  249396.199392857  293089.754507247


# # Add colname for first column
# col_names <- c("name", colnames(signal_corr_data)[-1])
# colnames(signal_corr_data) <- col_names

# Save as file and reload again
write.table(
    signal_corr_data,
    file=paste(output_path, "/signal_corr_data.tab", sep=""),
    sep="\t",
    row.names=TRUE,
    col.names=TRUE)

# Create sample list file
## sample batch class order
## 1 batch01_QC01 1  <NA>     1
## 2 batch01_QC02 1  <NA>     2

sample <- colnames(signal_corr_data[-1])
sample_idx <- match(sample, meta_all$file_name_pos)
batch <- meta_all[sample_idx, "block"]
class <- meta_all[sample_idx, "Regimen"]
order <- meta_all[sample_idx, "order"]
sampleListFile <- cbind(sample, batch, class, order)
colnames(sampleListFile) <- c("sample", "batch", "class", "order")
write.table(
    sampleListFile,
    file=paste(output_path, "/sampleListFile.tab", sep=""),
    sep="\t",
    row.names=TRUE,
    col.names=TRUE)
# sampleListFile <- read.table(file = "sampleListFile.tab", sep="\t")

# Do normalisation
para <- new("metaXpara")
pfile <- paste(output_path, "/signal_corr_data.tab", sep="")
sfile <- paste(output_path, "/sampleListFile.tab", sep="")
rawPeaks(para) <- read.delim(pfile, check.names = FALSE)
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
# para <- missingValueImpute(para)
# Why are there missing values appearing after QCRLSC normalisation???
para <- doQCRLSC(para, cpu=detectCores())
plotQCRLSC(para$metaXpara)

# Download data
para <- transformation(para$metaXpara, valueID="value")
non_norm_pos_peaklist <- getPeaksTable(para, valueID="value")
write.table(
    non_norm_pos_peaklist,
    file=paste(output_path, "/non_norm_pos_data.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)

para <- transformation(para$metaXpara, valueID = "valueNorm")
norm_pos_peaklist <- getPeaksTable(para, valueID="valueNorm")
# Format normalised negative peaklist data
# Use sample column as rownames
rownames(norm_pos_peaklist) <- norm_pos_peaklist$sample
# Remove columns
drop_cols <- c("sample", "class", "batch", "order")
norm_pos_peaklist <- norm_pos_peaklist[ , !(names(norm_pos_peaklist) %in% drop_cols)]
# Transpose data
norm_pos_peaklist <- t(norm_pos_peaklist)
# Order data based on the peak ID rownames
norm_pos_peaklist <- norm_pos_peaklist[naturalorder(rownames(norm_pos_peaklist)), ]
write.table(
    norm_pos_peaklist,
    file=paste(output_path, "norm_pos_data.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)

