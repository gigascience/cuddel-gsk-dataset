library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(specmine)
library(impute)
library(ggfortify)
library(metaX)

# Need to consider batch differences which might affect peaks in QCs due to the
# 4 analytical blocks are removed by XCMS during peak alignment.

######################################
# Read in metadata for data analysis #
######################################

# Read metadata file
meta_all <- read.csv("/home/peter/gsk/meta/meta_all.csv")
# Extract QCs
qc_meta_all <- subset(meta_all, type == "QC")
# Sort QCs by block then order
qc_meta_all <- qc_meta_all [ with(qc_meta_all, order(block, order)), ]
# Create file paths
neg_dir = "/home/peter/gsk/raw/esi_neg/netcdf"
get_neg_file_paths <- function(x, output) {
  file_name_neg <- x[3]
  type <- x[10]
  qc <- x[12]
  block <- x[13]
  if (block == "1") {
    file <- paste(neg_dir, "/block1neg/", file_name_neg, ".cdf", sep = "")
  }
  else if (block == "2") {
    file <- paste(neg_dir, "/block2neg/", file_name_neg, ".cdf", sep = "")
  }
  else if (block == "3") {
    file <- paste(neg_dir, "/block3neg/", file_name_neg, ".cdf", sep = "")
  }
  else if (block == "4") {
    file <- paste(neg_dir, "/block4neg/", file_name_neg, ".cdf", sep = "")
  }
}

neg_files <- apply(meta_all, 1, get_neg_file_paths, output = files)
# Sort negative QC filenames
neg_files <- naturalsort(neg_files)

########################
# Apply XCMS onto data #
########################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and will take time to complete!!
neg_xset <- xcmsSet(neg_files, step = 0.02, snthresh=3, mzdiff = 0.05)

# Match peaks representing same analyte across samples
grp_neg_xset <- group(neg_xset, bw = 10, mzwid = 0.05)

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

###############################################################
# Calculate percentage of missing values for each feature row #
###############################################################
percent_nas <-rep(0, nrow(neg_qc_peaklist))
calculate_percentage_na <- function(peak_intensities){
  for (i in 1:nrow(peak_intensities)) {
    feature_vector <- peak_intensities[i, ]
    na_row <- is.na(feature_vector)
    number_na <- sum(na_row)
    percent_na <- number_na/ncol(peak_intensities)*100
    percent_nas[i] <- percent_na
  }
  return(percent_nas)
}
percent_nas <- calculate_percentage_na(neg_qc_peaklist)

# Paste percent_nas onto QC data frame
neg_qc_peaklist <- cbind(neg_peaklist[1:12], percent_nas, neg_peaklist[neg_qc_fnames])

###############################################
# Remove rows with 40% or more missing values #
###############################################
remove_bad_feature_vectors <- function(peak_list, percentage_threshold = 40){
  rows_for_deleting <- vector('numeric')
  for (i in 1:nrow(peak_list)) {
    row <- peak_list[i, ]
    na_percentage <- row[, "percent_nas"]
    if (na_percentage > percentage_threshold) {
      #print(peak_list[i,])
      #print(peak_list[i,ncol(row)])
      rows_for_deleting <- c(rows_for_deleting, i)
    }
  }
  pretreated_peak_list <- peak_list[-rows_for_deleting, ]
  return(pretreated_peak_list)
}
del40_neg_qc_peaklist <- remove_bad_feature_vectors(neg_qc_peaklist)
del40_neg_qc_peaklist <- as.matrix(as.data.frame(lapply(del40_neg_qc_peaklist, as.numeric)))

###################################################################
# Impute means to replace remaining missing values using specmine #
###################################################################
# Create list object containing a data object that is the peaklist
# containing IDs as rownames and QC names are column names
specmine_qc_neg <- list(data = del40_neg_qc_peaklist[, 14:ncol(del40_neg_qc_peaklist)], type = "ms-spectra", description = "GSK pooled QC")

# Replace NAs with mean
impute_mean_qc_neg <- impute_nas_mean(specmine_qc_neg)

#############################################
# Do PCA plot to check effect of imputation #
#############################################
impute_mean_pca_data <- cbind(impute_mean_qc_neg$data)
rownames(impute_mean_pca_data) <- del40_neg_qc_peaklist[,"idx"]
impute_mean_pca_data <- t(impute_mean_pca_data)

# Create vector containing block information
sample_names <- colnames(del40_neg_qc_peaklist)
sample_names <- sample_names[14:length(sample_names)]
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
impute_mean_pca_data <- cbind(impute_mean_pca_data, block)
write.table(impute_mean_pca_data, file = "impute_mean_pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Now do PCA
impute_mean_pca_data <- read.table(file = "impute_mean_pca_data.csv", sep=",")
autoplot(prcomp(impute_mean_pca_data[,1:ncol(impute_mean_pca_data)-1]), data = impute_mean_pca_data, colour = 'block', main = 'PCA on negative QC data with imputed means for missing values')
ggsave("impute_mean_pca_data.png")

##################################################################################
# Impute k-nearest neighbours to replace remaining missing values using specmine #
##################################################################################

# Replace NAs with k-means clustering
impute_knn_qc_neg <- impute_nas_knn(specmine_qc_neg, k=10)

# Do PCA plot to check effects k-means values in the QCs from the 4
# analytical blocks
impute_knn_pca_data <- impute_knn_qc_neg$data
impute_knn_pca_data <- cbind(impute_knn_pca_data$data)
rownames(impute_knn_pca_data) <- del40_neg_qc_peaklist[,"idx"]
impute_knn_pca_data <- t(impute_knn_pca_data)
# Create vector containing block information
sample_names <- colnames(del40_neg_qc_peaklist)
sample_names <- sample_names[14:length(sample_names)]
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
impute_knn_pca_data <- cbind(impute_knn_pca_data, block)
write.table(impute_knn_pca_data, file = "impute_knn_pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Now do PCA
impute_knn_pca_data <- read.table(file = "impute_knn_pca_data.csv", sep=",")
autoplot(prcomp(impute_knn_pca_data[,1:ncol(impute_knn_pca_data)-1]), data = impute_knn_pca_data, colour = 'block', main = 'PCA on negative QC data with imputed K-means imputation for missing values')
ggsave("impute_knn_pca_data.png")


#####################################################
# Filter data containing k-means imputations by RSD #
#####################################################

# Prepare data for RSD analysis
impute_knn_qc_neg <- impute_nas_knn(specmine_qc_neg, k=10)
impute_knn_qc_neg <- impute_knn_qc_neg$data
impute_knn_qc_neg <- impute_knn_qc_neg$data
rownames(impute_knn_qc_neg) <- del40_neg_qc_peaklist[,"idx"]

# Process data by rows according to peaks across QC samples
# Peaks are removed if its RSD QC value was >20%; i.e. the analytical
# reproducibility of the peak was considered too high.
filter_by_rsd <- function(data, feature_percent_threshold = 20) {
  # Calculate mean for each feature vector row
  feature_vector_means <- rowMeans(data)
  # Calculate standard deviation by feature vector row
  feature_vector_standard_deviations <- apply(data, 1, sd)
  # Check feature_vector_means and feature_vector_standard_deviations lengths
  if (length(feature_vector_means) != length(feature_vector_standard_deviations)) {
  	print(paste0("feature_vector_means length is not equal to feature_vector_standard_deviations length"))
  }
  # Calculate RSD per row from above matrix and add as third column by dividing the standard deviation by the mean and then multiply the result by 100 to express it as a percentage.
  feature_vector_rsds <- apply(data[,1:84 ], 1, function(data) {
    sd <- sd(data)
	mean <- mean(data)
	rsd <- (sd/mean) * 100
	return(rsd)
  })

  # Column bind feature_vector_means and feature_vector_standard_deviations
  data <- cbind(data, feature_vector_means, feature_vector_standard_deviations,feature_vector_rsds)
  # List peak IDs from rownames which are above 20% RSD and need removing from
  # feature vector matrix
  rows_for_deleting <- which(data[,87] > 20)
  data <- data[-rows_for_deleting, ]
  return(data)
}

rsd_knn_qc_neg <- filter_by_rsd(impute_knn_qc_neg)

############################################################
# Do PCA to check effect of RSD filtering after imputation #
############################################################

# Create vector containing block information
sample_names <- colnames(rsd_knn_qc_neg)
sample_names <- sample_names[1:84]
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
rsd_knn_qc_neg <- cbind(rsd_knn_qc_neg, block)
write.table(rsd_knn_qc_neg, file = "rsd_knn_qc_neg.csv", sep =",", row.names = TRUE, col.names = TRUE)

rsd_knn_qc_neg_pca_data <- read.table(file = "rsd_knn_qc_neg.csv", sep=",")
autoplot(prcomp(rsd_knn_qc_neg_pca_data[,1:ncol(rsd_knn_qc_neg_pca_data)-1]), data = rsd_knn_qc_neg_pca_data, colour = 'block', main = 'PCA on RSD filtered, k-means imputed QC negative data')
ggsave("pca_rsd_knn_qc_neg.png")

###########################################
# Get sample file names which are not QCs #
###########################################
get_neg_sample_file_paths <- function(x, output) {
  file_name_neg <- x[3]
  type <- x[10]
  qc_meta_all <- x[12]
  block <- x[13]
  if (type == "Sample" & block == "1") {
    qc_file <- paste(neg_dir, "/block1neg/", file_name_neg, ".cdf", sep = "")
  }
  else if (type == "Sample" & block == "2") {
    qc_file <- paste(neg_dir, "/block2neg/", file_name_neg, ".cdf", sep = "")
  }
  else if (type == "Sample" & block == "3") {
    qc_file <- paste(neg_dir, "/block3neg/", file_name_neg, ".cdf", sep = "")
  }
  else if (type == "Sample" & block == "4") {
    qc_file <- paste(neg_dir, "/block4neg/", file_name_neg, ".cdf", sep = "")
  }
  else {
    qc_file <- "QC"
  }
}
neg_sample_files <- apply(meta_all, 1, get_neg_sample_file_paths, output = sample_files)
# Remove "QC" values in character vector
neg_sample_files <- neg_sample_files [! neg_sample_files %in% "QC"]
# Sort negative sample filenames
neg_sample_files <- naturalsort(neg_sample_files)

neg_sample_fnames <- substr(neg_sample_files, 46, 80)
neg_sample_fnames <- substr(neg_sample_fnames, 1, nchar(neg_sample_fnames)-4)

#####################################################
# Prepare negative QC and sample data for joint PCA #
#####################################################

# Combine QC and sample filenames
neg_fnames <- c(neg_qc_fnames, neg_sample_fnames)
# Create matrix containing QC and sample peaklists
neg_qc_sample_peaklist <- neg_peaklist[neg_fnames]

#################################################
# Do joint PCA to check negative peak list data #
#################################################
# Create vector containing QC and sample information
pca_meta_qc_sample <- integer(0)
pca_meta_qc_sample[1:length(neg_qc_fnames)] <- "QC"
pca_meta_qc_sample[length(neg_qc_fnames)+1:length(neg_sample_fnames)] <- "Sample"

# Create vector containing QC, sample and block information
pca_meta_qc_sample_block <- integer(0)
colnames_neg_qc_sample_peaklist <- colnames(neg_qc_sample_peaklist)
for (i in 1:length(colnames_neg_qc_sample_peaklist)) {
  if (colnames_neg_qc_sample_peaklist[i] %in% neg_qc_fnames && grepl("block1", colnames_neg_qc_sample_peaklist[i]) == 1) {
  	pca_meta_qc_sample_block[i] <- "QC_block1"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_qc_fnames && grepl("block2", colnames_neg_qc_sample_peaklist[i]) == 1) {
    pca_meta_qc_sample_block[i] <- "QC_block2"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_qc_fnames && grepl("block3", colnames_neg_qc_sample_peaklist[i]) == 1) {
      pca_meta_qc_sample_block[i] <- "QC_block3"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_qc_fnames && grepl("block4", colnames_neg_qc_sample_peaklist[i]) == 1) {
	pca_meta_qc_sample_block[i] <- "QC_block4"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_sample_fnames && grepl("block1", colnames_neg_qc_sample_peaklist[i]) == 1) {
	pca_meta_qc_sample_block[i] <- "Sample_block1"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_sample_fnames && grepl("block2", colnames_neg_qc_sample_peaklist[i]) == 1) {
  	pca_meta_qc_sample_block[i] <- "Sample_block2"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_sample_fnames && grepl("block3", colnames_neg_qc_sample_peaklist[i]) == 1) {
	pca_meta_qc_sample_block[i] <- "Sample_block3"
  }
  else if (colnames_neg_qc_sample_peaklist[i] %in% neg_sample_fnames && grepl("block4", colnames_neg_qc_sample_peaklist[i]) == 1) {
  	pca_meta_qc_sample_block[i] <- "Sample_block4"
  }
}

# Create vector containing block numbers for labelling shapes
pca_meta_qc_block <- integer(0)
colnames_neg_qc_sample_peaklist <- colnames(neg_qc_sample_peaklist)
for (i in 1:length(colnames_neg_qc_sample_peaklist)) {
  if (grepl("block1", colnames_neg_qc_sample_peaklist[i]) == 1) {
  	pca_meta_qc_block[i] <- "block1"
  }
  else if (grepl("block2", colnames_neg_qc_sample_peaklist[i]) == 1) {
    pca_meta_qc_block[i] <- "block2"
  }
  else if (grepl("block3", colnames_neg_qc_sample_peaklist[i]) == 1) {
      pca_meta_qc_block[i] <- "block3"
  }
  else if (grepl("block4", colnames_neg_qc_sample_peaklist[i]) == 1) {
	pca_meta_qc_block[i] <- "block4"
  }
}

# Transpose data
pca_data <- t(neg_qc_sample_peaklist)
# Add block information to PCA data
pca_data <- cbind(pca_data, pca_meta_qc_sample)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = 'pca_meta_qc_sample', main = 'PCA on unprocessed negative QC and sample data')
ggsave("unprocessed_neg_qc_sample_pca.png")

# Using different shapes and colours for labelling analytical blocks, QCs and samples
# Transpose data
pca_data <- t(neg_qc_sample_peaklist)
# Add block information to PCA data
pca_data <- cbind(pca_data, pca_meta_qc_sample, pca_meta_qc_block)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape = 'pca_meta_qc_block', colour = 'pca_meta_qc_sample', main = 'PCA on unprocessed negative QC and sample data')
ggsave("unprocessed_neg_qc_sample_pca.png")

#######################################################################
# Use metaX to perform data normalisation to correct within batch and #
# batch-to-batch variation                                            #
#######################################################################

# Use this data
ncol(rsd_knn_qc_neg)
[1] 87
length(colnames(rsd_knn_qc_neg))
[1] 87
# Remove not required columns
test = subset(rsd_knn_qc_neg, select = -c(feature_vector_means, feature_vector_standard_deviations, feature_vector_rsds))
# Move rownames into first column
test <- cbind(rownames(test), test)
# Check
head(test[,1:4])
# Add colname for first column
col_names <- c("name", colnames(test)[-1])
colnames(test) <- col_names
# Reformat row names
rownames(test) <- rep(1:nrow(test))
# Save as file and reload again
write.table(test, file = "test.tab", sep ="\t", row.names = TRUE, col.names = TRUE)
rawPeaks <- read.table(file = "test.tab", sep="\t")

# Create sample list file
## sample batch class order
## 1 batch01_QC01 1  <NA>     1
## 2 batch01_QC02 1  <NA>     2
sampleListFile <- cbind(rownames(meta_all), meta_all$block, meta_all$Regimen, meta_all$order)
colnames(sampleListFile) <- c("sample", "batch", "class", "order")
# Save as file and reload again
write.table(sampleListFile, file = "sampleListFile.tab", sep ="\t", row.names = TRUE, col.names = TRUE)
sampleListFile <- read.table(file = "sampleListFile.tab", sep="\t")

# Do normalisation
para <- new("metaXpara")
pfile <- "/home/peter/gsk/raw/esi_neg/netcdf/test.tab"
sfile <- "/home/peter/gsk/raw/esi_neg/netcdf/sampleListFile.tab"
rawPeaks(para) <- read.delim(pfile, check.names = FALSE)
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
# para <- missingValueImpute(para)
res <- doQCRLSC(para, cpu=1)
plotQCRLSC(res$metaXpara)

# Get peak data and do PCA
peaksData <- para@peaksData




