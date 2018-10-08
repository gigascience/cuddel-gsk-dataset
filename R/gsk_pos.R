# Title     : gsk_pos.R
# Objective : Processing of GSK positive mode data set
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

source("functions.R")

# Location of GSK data set
datadir = "/home/peter/"
# File path for positive files
pos_dir = paste(datadir, "gsk/raw/esi_pos/netcdf", sep="")
# Output directory
output_path <- paste(pos_dir, "/output", sep="")


######################################
# Read in metadata for data analysis #
######################################

# Read in and sort metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta_20180831.csv", sep=""))
meta <- meta[naturalorder(meta$file_name_pos), ]
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all_20180831.csv", sep=""))
meta_all <- meta_all[naturalorder(meta_all$file_name_pos), ]

# Get list of positive file paths
pos_file_paths <- getAllGSKFilePaths(mode="positive")

# Get list of positive files
pos_files <- meta[, "file_name_pos"]
pos_files <- as.character(pos_files)


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
# Get peaklist - contains 6200 peaks identified by XCMS
pos_peaklist <- getPeaklist(pos_xsa)
# Add numeric index to last column in pos_peaklist to keep track of peaks during
# data pre-treatment
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

pos_peaklist <- read.csv(paste(output_path, "/xcms_pos_peaklist.csv", sep=""))


#################################################
# Do PCA on unprocessed positive QC and samples #
#################################################

# Prepare pos_peaklist data
pca_data <- pos_peaklist[, pos_files]
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
    main='PCA scores plots for positive ion mode QC and plasma sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pos_qc_plasma_pca1.png", sep=""))

autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA scores plots for positive ion mode QC and plasma sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pos_qc_plasma_pca2.png", sep=""))

# There is an outlier by a block 4 plasma sample
# which is skewing the PCA result


##############################
# Do PCA on positive QC data #
##############################

# Prepare QC positive peaklist data
pos_qc_names <- getQCSampleNames(mode="positive")
# Subset QC data from positive peak list
pca_data <- pos_peaklist[, pos_qc_names]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
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
    main='PCA scores plots for positive ion mode QC sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pos_qc_pca.png", sep=""))


############################################################################
# Calculate percentage of missing values for each feature row for QCs only #
############################################################################

# Prepare QC positive data
pos_qc_peaklist <- pos_peaklist[, pos_qc_names]
# Calculate percentage of missing values by peak row
percent_nas <- rowMeans(is.na(pos_qc_peaklist))
# Add columns
idx <- rownames(pos_qc_peaklist)
pos_qc_peaklist <- cbind(idx, percent_nas, pos_qc_peaklist)
# Remove rows containing 40% or more missing values
na_filtered_pos_qc_peaklist <- pos_qc_peaklist[pos_qc_peaklist$percent_nas < .40, ]

# pos_qc_peaklist has 6200 rows. After filtering out 1747 rows with 40% na,
# na_filtered_pos_qc_peaklist contains 4453 rows


###############################################################
# Do PCA after removing features with 40% missing data values #
###############################################################

# Prepare positive peaklist data
unwanted_cols <- colnames(na_filtered_pos_qc_peaklist) %in% c("percent_nas", "idx")
pca_data <- na_filtered_pos_qc_peaklist[, !unwanted_cols]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
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
    main='PCA scores plots for positive ion mode QC samples after removing peaks with >40% NA',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/40_percent_filtered_pos_qc_pca.png", sep=""))


###################################################################
# Impute k-nearest neighbours to replace remaining missing values #
###################################################################

# Prepare data for knn analysis
unwanted_cols <- colnames(na_filtered_pos_qc_peaklist) %in% c("percent_nas", "idx")
impute_data <- na_filtered_pos_qc_peaklist[, !unwanted_cols]
# Transpose data
impute_data <- t(impute_data)  # 84 rows=samples, 4453 columns=peaks

# Identify which column peak features have missing values
na_peak_cols <- colnames(impute_data)[colSums(is.na(impute_data)) > 0]
# Perform K-means clustering using VIM - takes time!!
knn_pos_qc <- kNN(impute_data, variable=na_peak_cols)

# Output knn data because it's creation is compute-intensive
write.table(
    knn_pos_qc,
    file = paste(output_path, "/knn_pos_qc.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)

# Copy peak features from kNN result
knn_pos_qc_peaklist <- knn_pos_qc[, 1:ncol(impute_data)]
knn_pos_qc_peaklist <- as.matrix(knn_pos_qc_peaklist)
# Copy rownames and column names from impute_knn_data to rownames for knn
rownames(knn_pos_qc_peaklist) <- rownames(impute_data)
colnames(knn_pos_qc_peaklist) <- colnames(impute_data)
# Transpose
knn_pos_qc_peaklist <- t(knn_pos_qc_peaklist)


#####################################################
# Do PCA plot to check effect of imputed knn values #
#####################################################
pca_data <- knn_pos_qc_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
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
    main='PCA scores plot for positive QC data with imputed k-means for missing values',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/knn_pos_qc_pca.png", sep=""))


######################################
# Filter k-means imputed data by RSD #
######################################

# Prepare data for RSD analysis
rsd_qc_pos_peaklist <- knn_pos_qc_peaklist

# Calculate RSD
rowRSD <- apply(rsd_qc_pos_peaklist, 1, function(data) {
    sd <- sd(data)
    mean <- mean(data)
    rsd <- (sd/mean) * 100
    return(rsd)
})

# Find out which rows are over 20% RSD and remove them
bad_rows <- which(rowRSD > 20)  # there are 4430 peaks with >20% RSD
rsd_qc_pos_peaklist <- rsd_qc_pos_peaklist[-bad_rows, ]  # Only 23 peaks left!


###########################################
# Do PCA to check effect of RSD filtering #
###########################################

pca_data <- rsd_qc_pos_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
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
    main='PCA scores plot on RSD-filtered QC positive data',
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
# Only 1 peak remains after all rows containing missing values are removed.
# Therefore no point performing PCA.

# # Prepare metadata for labelling PCA graph
# samples <- colnames(pca_data)
# block_meta <- getBlockMetadata(samples, mode="positive")
# # Sample type information
# sample_type_meta <- getQCSampleMetadata(samples, mode="positive")
# # Convert block_meta character vector to matrix
# pca_meta <- cbind(block_meta, sample_type_meta)
#
# # Transpose data
# pca_data <- t(pca_data)
#
# # Plot PCA results
# autoplot(
#     prcomp(pca_data),
#     data=pca_meta,
#     shape='sample_type_meta',
#     colour='block_meta',
#     main='PCA on filtered positive QCs with positive plasma samples',
#     frame=TRUE,
#     frame.type='norm')
# ggsave(paste(output_path, "/filtered_pos_qc_plasma_samples_pca.png", sep=""))
#
# # Plot PCA results
# autoplot(
#     prcomp(pca_data),
#     data=pca_meta,
#     shape='block_meta',
#     colour='sample_type_meta',
#     main='PCA on pretreated positive QCs with positive plasma samples',
#     frame=TRUE,
#     frame.type='norm')
# ggsave(paste(output_path, "/filtered_pos_qc_plasma_samples_pca_shape_block_colour_sample.png", sep=""))


###########################################################
# Use k-means to impute missing values for plasma samples #
###########################################################

# Prepare data for imputation
impute_data <- filtered_qc_plasma_peaklist
impute_data <- t(impute_data)

# Identify column peak features with missing values
na_cols <- colnames(impute_data)[colSums(is.na(impute_data)) > 0]
# Perform K-means clustering using VIM
knn_pos_qc_plasma <- kNN(impute_data, variable=na_cols)

# Copy peak features from kNN result
knn_qc_plasma_peaklist <- knn_pos_qc_plasma[, 1:ncol(impute_data)]
knn_qc_plasma_peaklist <- as.matrix(knn_qc_plasma_peaklist)
# Copy rownames and column names
rownames(knn_qc_plasma_peaklist) <- rownames(impute_data)
colnames(knn_qc_plasma_peaklist) <- colnames(impute_data)

# Transpose
knn_qc_plasma_peaklist <- t(knn_qc_plasma_peaklist)


#########################################################################
# Do PCA on combined rsd_qc_pos_peaklist and pretreated_sample_peaklist #
#########################################################################

pca_data <- knn_qc_plasma_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
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
    main='PCA scores plot for pretreated positive QCs with knn positive samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pretreated_qc_knn_pos_samples_pca1.png", sep=""))

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='sample_type_meta',
    colour='block_meta',
    main='PCA scores plot for pretreated positive QCs with knn positive samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pretreated_qc_knn_pos_samples_pca2.png", sep=""))


###########################################################################
# Use metaX data normalisation to correct within batch and batch-to-batch #
# variation in positive QCs and samples                                   #
###########################################################################

# Need to consider batch differences affecting peaks in QCs due to the 4
# analytical blocks are removed by XCMS during peak alignment.

# Prepare data
signal_corr_data <- knn_qc_plasma_peaklist
signal_corr_data <- cbind(rownames(signal_corr_data), signal_corr_data)
colnames(signal_corr_data)[1] <- "name"
signal_corr_data <- data.frame(signal_corr_data)

# Check
head(signal_corr_data[, 1:4])
##     name GSK_pos_block1_09 GSK_pos_block1_10 GSK_pos_block1_16
## 210  210  203705.219952001  220888.421936538  179189.455627451
## 257  257   106059.74673913      85431.038875  102438.816159574
## 439  439  177268.260333333   177413.66791746  177268.260333333
## 478  478  319688.213197872        238508.938  319688.213197872
## 727  727  85638.9250857143  85638.9250857143  84019.8217129621
## 815  815  107910.199359322     130026.891285     130026.891285


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
# Format normalised positive peaklist data
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
    file=paste(output_path, "/norm_pos_data.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)


#####################################
# Do PCA to check signal correction #
#####################################

# Prepare data
pca_data <- norm_pos_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
sample_type_meta <- getQCSampleMetadata(samples, mode="positive")
regimen_meta <- getRegimenMetadata(samples, mode="positive")
# Convert block_meta character vector to matrix
pca_meta <- cbind(block_meta, sample_type_meta, regimen_meta)

# Transpose data
pca_data <- t(pca_data)

# Highlight sample type and batch in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_type_meta',
    main='PCA scores plot for normalised positive QC and plasma sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/norm_pos_qc_sample_pca1.png", sep=""))

# Highlight regimens in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    colour='regimen_meta',
    main='PCA scores plot for normalised positive QC and plasma sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/norm_pos_qc_sample_pca2.png", sep=""))

# Highlight blocks in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    colour='block_meta',
    shape='sample_type_meta',
    main='PCA scores plot for normalised positive QC and plasma sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/norm_pos_qc_sample_pca3.png", sep=""))


############################################################################
# Normalise plasma data to QCs by division of the median feature intensity #
# responses measured for QC samples with the intensity response of each    #
# feature for a plasma sample                                              #
############################################################################

# Get all positive QC sample names
pos_qc_names <- getQCSampleNames(mode="positive")
# Get positive QC data
pos_norm_qc_data <- norm_pos_peaklist[, pos_qc_names]
# Calculate median peak intensity
pos_qc_medians <- apply(pos_norm_qc_data, 1, median, na.rm=TRUE)
# Get all positive plasma sample names
pos_plasma_names <- getPlasmaSampleNames(mode="positive")
# Get positive QC data
pos_norm_plasma_peaklist <- norm_pos_peaklist[, pos_plasma_names]
# Divide all plasma peak intensity values with the relevant median QC peak intensity
median_norm_pos_plasma_peaklist <- pos_norm_plasma_peaklist/pos_qc_medians


###########################################
# Do PCA to check QC median normalisation #
###########################################

# Prepare data containing both QC and plasma sample data
pca_data <- cbind(pos_norm_qc_data, median_norm_pos_plasma_peaklist)

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
sample_type_meta <- getQCSampleMetadata(samples, mode="positive")
regimen_meta <- getRegimenMetadata(samples, mode="positive")
pca_meta <- cbind(block_meta, sample_type_meta, regimen_meta)

# Transpose data
pca_data <- t(pca_data)

# Highlight sample type and batch in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    colour='regimen_meta',
    main='PCA scores plot for QC-median normalised positive plasma data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/QC_median_norm_pos_QC_plasma_pca.png", sep=""))


##############################################################################
# Normalise data to time point 0 to partially compensate for metabolic       #
# differences among the volunteers and to identify the changes in metabolome #
# that are relative to the baseline                                          #
##############################################################################

# Get all TP 0 data
tp0_pos_plasma_names <- getTP0PlasmaSampleNames(mode="positive")
tp0_pos_plasma_peaklist <- norm_pos_peaklist[, tp0_pos_plasma_names]
# Calculate average value for each TP0 peak
tp0_means <- rowMeans(tp0_pos_plasma_peaklist)
# Divide all data with tp0 means
tp0_norm_plasma_peaklist <- median_norm_pos_plasma_peaklist/tp0_means

# Download data as CSV file for combining with negative data
write.table(tp0_norm_plasma_peaklist, file="pos_tp0_norm_plasma_peaklist.csv", sep=",", row.names = TRUE, col.names = TRUE)


#####################################
# Do PCA to check TP0 normalisation #
#####################################

# Prepare data
pca_data <- cbind(tp0_norm_plasma_peaklist)

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples, mode="positive")
sample_type_meta <- getQCSampleMetadata(samples, mode="positive")
regimen_meta <- getRegimenMetadata(samples, mode="positive")
pca_meta <- cbind(block_meta, sample_type_meta, regimen_meta)

# Transpose data
pca_data <- t(pca_data)

# Highlight sample type and batch in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    colour='regimen_meta',
    main='PCA scores plot for TP0 and QC-median normalised positive plasma sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/TP0_norm_neg_peaklist_pca.png", sep=""))
