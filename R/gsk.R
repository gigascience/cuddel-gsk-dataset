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
# File path for negative files
neg_dir = paste(datadir, "gsk/raw/esi_neg/netcdf", sep ="")

# Output path
output_path <- paste(neg_dir, "/output", sep="")

######################################
# Read in metadata for data analysis #
######################################

# Read in and sort metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta.csv", sep =""))
meta <- meta[naturalorder(meta$file_name_neg),]
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all.csv", sep =""))
meta_all <- meta_all[naturalorder(meta_all$file_name_neg),]

# Get list of negative file paths
neg_file_paths <- getAllNegGSKFilePaths()

# Get list of negative files
neg_files <- meta_all[,"file_name_neg"]
neg_files <- as.character(neg_files)

######################################
# Apply XCMS onto QC and sample data #
######################################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and takes time to complete!!
neg_xset <- xcmsSet(neg_file_paths, step=0.02, snthresh=3, mzdiff=0.05)

# Match peaks representing same analyte across samples
grp_neg_xset <- group(neg_xset, bw=10, mzwid=0.05)

# Create CAMERA object
neg_xsa <- xsAnnotate(grp_neg_xset)
# Get peaklist
neg_peaklist <- getPeaklist(neg_xsa)
# Add numeric index to negative peaklist to keep track of peaks during data pre-treatment
neg_peaklist$idx <- seq.int(nrow(neg_peaklist))
# Move index to left hand side of data frame
neg_peaklist <- neg_peaklist[, c(ncol(neg_peaklist), 1:(ncol(neg_peaklist)-1))]

# Output neg_peaklist data
write.table(
    neg_peaklist,
    file = paste(output_path, "/xcms_neg_peaklist.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)

#############################
# Load saved peak list data #
#############################

xcms_neg_peaks = read.csv(paste(output_path, "/xcms_neg_peaklist.csv", sep=""))


#################################################
# Do PCA on unprocessed negative QC and samples #
#################################################

# Prepare neg_peaklist data
pca_data <- xcms_neg_peaks[, neg_files]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
sample_meta <- getQCSampleMetadata(samples)
pca_meta <- cbind(block_meta, sample_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_meta',
    main='PCA on unprocessed negative QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_neg_qc_sample_pca.png", sep=""))

autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on unprocessed negative QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_neg_qc_sample_pca.png", sep=""))


##############################
# Do PCA on negative QC data #
##############################

# Prepare QC negative peaklist data
meta_all_qc_rows <- meta_all[, "type"]=="QC"
neg_qc_names <- meta_all[meta_all_qc_rows, "file_name_neg"]
neg_qc_names <- as.character(neg_qc_names)
# Subset QC data from negative peak list
pca_data <- xcms_neg_peaks[, neg_qc_names]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on unprocessed negative QCs',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/unprocessed_neg_qc_pca.png", sep=""))


############################################################################
# Calculate percentage of missing values for each feature row for QCs only #
############################################################################

# Prepare QC negative data
qc_neg_peaklist <- xcms_neg_peaks[, neg_qc_names]
# Calculate percentage of missing values by peak row
percent_nas <- rowMeans(is.na(qc_neg_peaklist))
# Add columns
idx <- rownames(qc_neg_peaklist)
qc_neg_peaklist <- cbind(idx, percent_nas, qc_neg_peaklist)
# Remove rows containing 40% or more missing values
na_filtered_qc_neg_peaklist <- qc_neg_peaklist[qc_neg_peaklist$percent_nas < .40, ]

# qc_neg_peaklist has 3134 rows. After filtering rows with 40% na, na_filtered_qc_neg_peaklist contains 1704 rows


###############################################################
# Do PCA after removing features with 40% missing data values #
###############################################################

# Prepare negative peaklist data
nr <- colnames(na_filtered_qc_neg_peaklist) %in% c("percent_nas", "idx")
pca_data <- na_filtered_qc_neg_peaklist[, !nr]
# Remove row peaks containing missing values
pca_data <- pca_data[rowSums(is.na(pca_data))==0, ]

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on negative QCs after peaks missing 40% values removed',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/40_percent_filtered_neg_qc_pca.png", sep=""))


###################################################################
# Impute k-nearest neighbours to replace remaining missing values #
###################################################################

# Prepare data for knn analysis
nr <- colnames(na_filtered_qc_neg_peaklist) %in% c("percent_nas", "idx")
impute_data <- na_filtered_qc_neg_peaklist[, !nr]
# Transpose data
impute_data <- t(impute_data)

# Identify which column peak features have missing values
na_peak_cols <- colnames(impute_data)[colSums(is.na(impute_data)) > 0]
# Perform K-means clustering using VIM - takes time!!
knn <- kNN(impute_data, variable=na_peak_cols)
# Copy peak features from kNN result
knn_qc_neg_peaklist <- knn[, 1:ncol(impute_data)]
knn_qc_neg_peaklist <- as.matrix(knn_qc_neg_peaklist)
# Copy rownames and column names from impute_knn_data to rownames for knn
rownames(knn_qc_neg_peaklist) <- rownames(impute_data)
colnames(knn_qc_neg_peaklist) <- colnames(impute_data)
# Transpose
knn_qc_neg_peaklist <- t(knn_qc_neg_peaklist)


#####################################################
# Do PCA plot to check effect of imputed knn values #
#####################################################
pca_data <- knn_qc_neg_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on negative QC data with imputed k-means for missing values',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/impute_knn_pca_data.png", sep=""))


######################################
# Filter k-means imputed data by RSD #
######################################

# Prepare data for RSD analysis
rsd_qc_neg_peaklist <- knn_qc_neg_peaklist

# Calculate RSD
rowRSD <- apply(rsd_qc_neg_peaklist, 1, function(data) {
  sd <- sd(data)
  mean <- mean(data)
  rsd <- (sd/mean) * 100
  return(rsd)
})

# Find out which rows are over 20% RSD and remove them
bad_rows <- which(rowRSD > 20)
rsd_qc_neg_peaklist <- rsd_qc_neg_peaklist[-bad_rows, ]


###########################################
# Do PCA to check effect of RSD filtering #
###########################################

pca_data <- rsd_qc_neg_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='block_meta',
    main='PCA on RSD filtered QC negative data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/rsd_neg_qc_pca_data.png", sep=""))


##########################################################################
# Combine feature vectors from plasma samples with filtered negative QCs #
##########################################################################

# Get peak indices of pretreated QC feature vectors
pretreated_qc_idx <- rownames(rsd_qc_neg_peaklist)
# Get plasma peaklist from neg_peaklist based on pretreated QC indices
filtered_neg_peaklist <- neg_peaklist[rownames(neg_peaklist) %in% pretreated_qc_idx, ]
neg_plasma_sample_names <- getNegPlasmaSampleNames()
filtered_plasma_peaklist <- filtered_neg_peaklist[, colnames(filtered_neg_peaklist) %in% neg_plasma_sample_names]
# Combine QC and plasma samples from negative peaklist
filtered_qc_plasma_peaklist <- cbind(rsd_qc_neg_peaklist, filtered_plasma_peaklist)


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
sample_type_meta <- getQCSampleMetadata(samples)
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta, sample_type_meta)

# Transpose data
pca_data <- t(pca_data)

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='sample_type_meta',
    colour='block_meta',
    main='PCA on filtered negative QCs with negative plasma samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/filtered_neg_qc_plasma_samples.png", sep=""))

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_type_meta',
    main='PCA on pretreated negative QCs with negative plasma samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/filtered_neg_qc_plasma_samples.png", sep=""))


###########################################################
# Use k-means to impute missing values for plasma samples #
###########################################################

# Prepare data for imputation
impute_data <- filtered_qc_plasma_peaklist
impute_data <- t(impute_data)

# Identify  column peak features with missing values
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
sample_type_meta <- getQCSampleMetadata(samples)
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
    main='PCA on pretreated negative QCs with knn negative samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pretreated_qc_knn_neg_samples_pca.png", sep=""))

# Plot PCA results
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_type_meta',
    main='PCA on pretreated negative QCs with knn negative samples',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/pretreated_qc_knn_neg_samples_pca.png", sep=""))


############################################################
# Use metaX data normalisation to correct within batch and #
# batch-to-batch variation in negative QCs and samples     #
############################################################

# Need to consider batch differences affecting peaks in QCs due to the
# 4 analytical blocks are removed by XCMS during peak alignment.

# Data for metaX signal correction needs to look like this:
# head(signal_correction_data[, 1:4])
##        name batch01_QC01 batch01_QC02 batch01_QC03
## 1  78.02055      14023.0      13071.0      15270.0
## 2 452.00345      22455.0      10737.0      27397.0
## 3 138.96337       6635.4       8062.3       6294.6
## 4  73.53838      26493.0      26141.0      25944.0
## 5 385.12885      57625.0      56964.0      59045.0
## 6 237.02815     105490.0      90166.0      92315.0

# Prepare data
signal_corr_data <- knn_qc_plasma_peaklist
signal_corr_data <- cbind(rownames(signal_corr_data), signal_corr_data)
colnames(signal_corr_data)[1] <- "name"
signal_corr_data <- data.frame(signal_corr_data)
# Check
head(signal_corr_data[, 1:4])
## name GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_16r
## 2     2   2214161.06684043   2213167.30024468   2200922.27885106
## 11   11   80404.2911244445   90185.9247093023   80860.0580355555
## 12   12   1052155.88762273   1013478.89673721   968455.803136363
## 46   46   95204.1498043478   95296.5965466666   96498.8371822222
## 47   47    37436017.503675   36930624.0750681   37326859.8432458
## 51   51         47954.2056   43746.2516363636   53033.2722418605

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
sample_idx <- match(sample, meta_all$file_name_neg)
batch <- meta_all[sample_idx, "block"]
class <- meta_all[sample_idx, "Regimen"]
order <- meta_all[sample_idx, "order"]
sampleListFile <- cbind(sample, batch, class, order)
colnames(sampleListFile) <- c("sample", "batch", "class", "order")
write.table(
    sampleListFile,
    file=paste(output_path, "/sampleListFile.tab", sep=""),
    sep ="\t",
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
para <- doQCRLSC(para, cpu=detectCores())
plotQCRLSC(para$metaXpara)

# Download data
para <- transformation(para$metaXpara, valueID="value")
non_norm_neg_peaklist <- getPeaksTable(para, valueID="value")
write.table(
    non_norm_neg_data,
    file=paste(output_path, "/non_norm_neg_data.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)

para <- transformation(para$metaXpara, valueID = "valueNorm")
norm_neg_peaklist <- getPeaksTable(para, valueID="valueNorm")
# Format normalised negative peaklist data
# Use sample column as rownames
rownames(norm_neg_peaklist) <- norm_neg_peaklist$sample
# Remove columns
drop_cols <- c("sample", "class", "batch", "order")
norm_neg_peaklist <- norm_neg_peaklist[ , !(names(norm_neg_peaklist) %in% drop_cols)]
# Transpose data
norm_neg_peaklist <- t(norm_neg_peaklist)
# Order data based on the peak ID rownames
norm_neg_peaklist <- norm_neg_peaklist[naturalorder(rownames(norm_neg_peaklist)), ]
write.table(
    norm_neg_data,
    file=paste(output_path, "norm_neg_data.csv", sep=""),
    sep=",",
    row.names=TRUE,
    col.names=TRUE)


#####################################
# Do PCA to check signal correction #
#####################################

# Prepare data
pca_data <- norm_neg_peaklist

# Prepare metadata for labelling PCA graph
samples <- colnames(pca_data)
block_meta <- getBlockMetadata(samples)
# Sample type information
sample_type_meta <- getQCSampleMetadata(samples)
regimen_meta <- getRegimenMetadata(samples)
# Convert to block_meta character vector to matrix
pca_meta <- cbind(block_meta, sample_type_meta, regimen_meta)

# Transpose data
pca_data <- t(pca_data)

# Highlight sample type and batch in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    shape='block_meta',
    colour='sample_type_meta',
    main='PCA on normalised negative QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/norm_neg_qc_sample_pca_coloured_sample_type.png", sep=""))

# Highlight regimens in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    colour='regimen_meta',
    main='PCA on normalised negative QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/norm_neg_qc_sample_pca_coloured_regimens.png", sep=""))

# Highlight blocks in PCA plot
autoplot(
    prcomp(pca_data),
    data=pca_meta,
    colour='block_meta',
    shape='sample_type_meta',
    main='PCA on normalised negative QC and sample data',
    frame=TRUE,
    frame.type='norm')
ggsave(paste(output_path, "/norm_neg_qc_sample_pca_coloured_blocks.png", sep=""))


#######################################################
# Normalise data to QCs by division of the median of  #
# feature intensity responses measured for QC samples #
# with the intensity response of each feature for a   #
# normal sample                                       #
#######################################################

# Get all negative QC sample names
neg_qc_names <- getNegativeQCSampleNames()
# Get negative QC data
neg_norm_qc_data <- norm_neg_peaklist[, qc_names]
# Calculate median peak intensity
neg_qc_medians <- apply(neg_norm_qc_data, 1, median, na.rm = TRUE)
# Get all negative plasma sample names
neg_plasma_names <- getNegativePlasmaSampleNames()
# Get negative QC data
neg_norm_plasma_data <- norm_neg_peaklist[, neg_plasma_names]
# Divide all plasma peak intensity values with the relevant median QC peak intensity



##################################
# Normalise data to time point 0 #
##################################

# Get all TP 0 data
tp0_sample_names <- meta_all[meta_all[, "Timepoint..PIMS."] == 0, "file_name_neg"]
tp0_sample_names <- as.character(tp0_sample_names)
tp0_sample_names <- na.omit(tp0_sample_names)
tp0_data <- norm_neg_peaklist[, tp0_sample_names]

# Calculate average value for each TP0 peak
tp0_means <- rowMeans(tp0_data)

# Divide all data with tp0 mean
tp0_qc_sample_neg <- cbind(qcnorm_qc_sample_neg, tp0_means)
tp0_qc_sample_neg <- as.matrix(tp0_qc_sample_neg)
counter <- 0
ans <- apply(tp0_qc_sample_neg, 1, function(x) {
  counter <<- counter + 1
  # print(x)
  # print(paste("Length of x:", length(x)))
  tp0_average <- as.numeric(x[length(x)])
  # print("##### End of matrix row #####")
  # for(i in 1:length(x)-1) {
  # print(x[length(x)])
  #   print(x[i] / x[length(x)])
  # }
  stuff <- lapply(x, function (y) {
    # print(paste("y:", y))
    # qc_average <- as.numeric(y[length(y)])
    # print(paste("qc_average: ", qc_average))
    new_value <- y / tp0_average
    # print(paste("new value: ", new_value))
    # print("#####")
    return(new_value)
  })
  # print(paste("stuff: ", stuff))
  # print(stuff)
  stuff <- as.numeric(stuff)
  # print(paste("length of stuff: ", length(stuff)))
  # print(class(stuff))
  # print(paste("counter:", counter))
  tp0_qc_sample_neg[counter,] <<- stuff
})

# Clean up
tp0_qc_sample_neg <- tp0_qc_sample_neg[, -ncol(tp0_qc_sample_neg)]

# Do PCA
pca_data <- tp0_qc_sample_neg

# Transpose data
pca_data <- t(pca_data)
# Add regimen information to PCA data
pca_data <- cbind(pca_data, regimen)

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = 'regimen', main = 'PCA on TP0-normalised negative QC and sample data', frame.type = 'norm')
ggsave("pca_qc_TP0_normalised.png")