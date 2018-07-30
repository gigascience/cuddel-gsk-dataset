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

# Let R know where GSK data set is located
datadir = "/home/peter/"

######################################
# Read in metadata for data analysis #
######################################

# Read metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta.csv", sep =""))
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all.csv", sep =""))
# Sort by file name neg
meta <- meta[naturalorder(meta$file_name_neg),]
meta_all <- meta_all[naturalorder(meta_all$file_name_neg),]

# Create file paths for negative files
neg_dir = paste(datadir, "gsk/raw/esi_neg/netcdf", sep ="")
file_name_neg <- meta_all[,"file_name_neg"]
file_name_neg <- as.character(file_name_neg)
neg_file_paths <- character()
block <- meta_all[,"block"]
for (i in 1:length(file_name_neg)) {
  if (block[i] == 1) {
    neg_file_path <- paste(neg_dir, "/block1neg/", file_name_neg[i], ".cdf", sep = "")
    neg_file_paths[i] <- neg_file_path
  }
  else if (block[i] == 2) {
    neg_file_path <- paste(neg_dir, "/block2neg/", file_name_neg[i], ".cdf", sep = "")
    neg_file_paths[i] <- neg_file_path
  }
  else if (block[i] == 3) {
    neg_file_path <- paste(neg_dir, "/block3neg/", file_name_neg[i], ".cdf", sep = "")
    neg_file_paths[i] <- neg_file_path
  }
  else if (block[i] == 4) {
    neg_file_path <- paste(neg_dir, "/block4neg/", file_name_neg[i], ".cdf", sep = "")
    neg_file_paths[i] <- neg_file_path
  }
}

######################################
# Apply XCMS onto QC and sample data #
######################################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and takes time to complete!!
neg_xset <- xcmsSet(neg_file_paths, step = 0.02, snthresh=3, mzdiff = 0.05)

# Match peaks representing same analyte across samples
grp_neg_xset <- group(neg_xset, bw = 10, mzwid = 0.05)

# Create CAMERA object
neg_xsa <- xsAnnotate(grp_neg_xset)
# Get peaklist
neg_peaklist <- getPeaklist(neg_xsa)
# Add numeric index to negative peaklist to keep track of peaks during data pre-treatment
neg_peaklist$idx <- seq.int(nrow(neg_peaklist))
# Move index to left hand side of data frame
neg_peaklist <- neg_peaklist[, c(ncol(neg_peaklist), 1:(ncol(neg_peaklist)-1))]

# Output neg_peaklist data
write.table(neg_peaklist, file = "neg_peaklist.csv", sep =",", row.names = TRUE, col.names = TRUE)

#############################
# Load saved peak list data #
#############################

xcms_neg_peaks = read.csv("xcms_neg_peaks.csv")


#################################################
# Do PCA on unprocessed negative QC and samples #
#################################################

# Prepare neg_peaklist data
pca_data <- xcms_neg_peaks[, file_name_neg]

# Prepare block information for labelling data points
block <- integer(0)
qc_sample_names <- colnames(pca_data)
for (i in 1:length(qc_sample_names)) {
  if (grepl("block1", qc_sample_names[i]) == 1) {
    block[i] <- "block1"
  }
  else if (grepl("block2", qc_sample_names[i]) == 1) {
    block[i] <- "block2"
  }
  else if (grepl("block3", qc_sample_names[i]) == 1) {
    block[i] <- "block3"
  }
  else if (grepl("block4", qc_sample_names[i]) == 1) {
    block[i] <- "block4"
  }
}

# Prepare QC and sample information for labelling data points
pca_meta_qc_sample <- integer(0)
for (i in 1:length(qc_sample_names)) {
  if (meta_all[which(meta_all[,"file_name_neg"]==qc_sample_names[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

# Transpose data
pca_data <- t(pca_data)
# Add block information, QC and sample names to PCA data
pca_data <- cbind(pca_data, block, pca_meta_qc_sample)
# PCA cannot be performed on data with missing values
# Therefore remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape= "block", colour = 'pca_meta_qc_sample', main = 'PCA on unprocessed negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("unprocessed_neg_qc_sample_pca.png")

autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape= "block", colour = 'block', main = 'PCA on unprocessed negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("unprocessed_neg_qc_sample_pca.png")


###################################################
# Extract negative QC data from negative peaklist #
###################################################

# Get QC sample names
meta_all_qc_rows <- meta_all[, "type"] == "QC"
neg_qc_names <- meta_all[meta_all_qc_rows, "file_name_neg"]
neg_qc_names <- as.character(neg_qc_names)
# Subset QC data from negative peak list
qc_neg_peaklist <- xcms_neg_peaks[, neg_qc_names]


##############################
# Do PCA on negative QC data #
##############################

# Prepare neg_peaklist data
pca_data <- qc_neg_peaklist

# Prepare block information for labelling data points
block <- integer(0)
qc_names <- colnames(pca_data)
for (i in 1:length(qc_names)) {
  if (grepl("block1", qc_names[i]) == 1) {
    block[i] <- "block1"
  }
  else if (grepl("block2", qc_names[i]) == 1) {
    block[i] <- "block2"
  }
  else if (grepl("block3", qc_names[i]) == 1) {
    block[i] <- "block3"
  }
  else if (grepl("block4", qc_names[i]) == 1) {
    block[i] <- "block4"
  }
}

# Transpose data
pca_data <- t(pca_data)
# Add block information, QC and sample names to PCA data
pca_data <- cbind(pca_data, block)
# PCA cannot be performed on data with missing values
# Therefore remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape = 'block', colour= "block", main = 'PCA on unprocessed negative QCs', frame = TRUE, frame.type = 'norm')
ggsave("unprocessed_neg_qc_pca.png")


############################################################################
# Calculate percentage of missing values for each feature row for QCs only #
############################################################################

percent_na_qc_neg_peaklist <- qc_neg_peaklist

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

percent_nas <- rep(0, nrow(percent_na_qc_neg_peaklist))
percent_nas <- calculate_percentage_na(percent_na_qc_neg_peaklist)

# Paste percent_nas into negative peak list data frame
percent_na_qc_neg_peaklist <- cbind(percent_na_qc_neg_peaklist, percent_nas)
# Copy peak id rownames into last column
percent_na_qc_neg_peaklist <- cbind(percent_na_qc_neg_peaklist, rownames(percent_na_qc_neg_peaklist))
colnames(percent_na_qc_neg_peaklist)[ncol(percent_na_qc_neg_peaklist)] <- "idx"

###############################################
# Remove rows with 40% or more missing values #
###############################################

over_40_percent_rows <- as.numeric(percent_na_qc_neg_peaklist$percent_nas) < 40
del40_qc_neg_peaklist <- percent_na_qc_neg_peaklist[over_40_percent_rows, ]

# qc_neg_peaklist contains 3134 rows. After filtering rows with 40% na = 1704 rows


###############################################################
# Do PCA after removing features with 40% missing data values #
###############################################################

# Not required columns for PCA analysis
nr <- c("percent_nas", "idx")
nr <- colnames(del40_qc_neg_peaklist) %in% c("percent_nas", "idx")
# Prepare neg_peaklist data
pca_data <- del40_qc_neg_peaklist[, !nr]

# Prepare block information for labelling data points
block <- integer(0)
qc_names <- colnames(pca_data)
for (i in 1:length(qc_names)) {
  if (grepl("block1", qc_names[i]) == 1) {
    block[i] <- "block1"
  }
  else if (grepl("block2", qc_names[i]) == 1) {
    block[i] <- "block2"
  }
  else if (grepl("block3", qc_names[i]) == 1) {
    block[i] <- "block3"
  }
  else if (grepl("block4", qc_names[i]) == 1) {
    block[i] <- "block4"
  }
}

# Transpose data
pca_data <- t(pca_data)
# Add block information to PCA data
pca_data <- cbind(pca_data, block)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, shape = 'block', colour= "block", main = 'PCA on negative QCs after peaks missing 40% values removed', frame = TRUE, frame.type = 'norm')
ggsave("40_percent_filtered_neg_qc_pca.png")


####################################################################
# Impute means to replace remaining missing values in negative QCs #
####################################################################

# Prepare  data
nr <- c("percent_nas", "idx")
nr <- colnames(del40_qc_neg_peaklist) %in% c("percent_nas", "idx")
impute_mean_qc_neg_peaklist <- del40_qc_neg_peaklist[,!nr]
na_rows <- which(is.na(impute_mean_qc_neg_peaklist), arr.ind=TRUE)
impute_mean_qc_neg_peaklist[na_rows] <- rowMeans(impute_mean_qc_neg_peaklist, na.rm=TRUE)[na_rows[,1]]


################################################
# Do PCA plot to check effect of imputed means #
################################################
pca_data <- impute_mean_qc_neg_peaklist
pca_data <- t(pca_data)

# Create vector containing block information
sample_names <- rownames(pca_data)
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

pca_data <- cbind(pca_data, block)
write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Now do PCA
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[1:ncol(pca_data)-1]), data = pca_data, colour= "block", main = 'PCA on negative QC data with imputed means for missing values')
ggsave("impute_mean_pca_data.png")


###################################################################
# Impute k-nearest neighbours to replace remaining missing values #
###################################################################

# Replace NAs with k-means clustering
impute_knn_qc_neg_peaklist <- del40_qc_neg_peaklist
# Clean up columns
nr <- c("percent_nas", "idx")
nr <- colnames(del40_qc_neg_peaklist) %in% c("percent_nas", "idx")
impute_knn_qc_neg_peaklist <- impute_knn_qc_neg_peaklist[, !nr]
# Transpose data
impute_knn_qc_neg_peaklist <- t(impute_knn_qc_neg_peaklist)
# Identify which column peak features have missing values
na_cols <- colnames(impute_knn_qc_neg_peaklist)[colSums(is.na(impute_knn_qc_neg_peaklist)) > 0]
# Perform K-means clustering using VIM
knn <- kNN(impute_knn_qc_neg_peaklist, variable = na_cols)
# Copy peak features from kNN result
impute_knn_qc_neg_peaklist2 <- knn[, 1:1704]
impute_knn_qc_neg_peaklist2 <- as.matrix(impute_knn_qc_neg_peaklist2)
# Copy rownames from impute_knn_data to rownames for knn
rownames(impute_knn_qc_neg_peaklist2) <- rownames(impute_knn_qc_neg_peaklist)
# Copy colnames from impute_knn_data to colnames for knn
colnames(impute_knn_qc_neg_peaklist2) <- colnames(impute_knn_qc_neg_peaklist)


#####################################################
# Do PCA plot to check effect of imputed knn values #
#####################################################
pca_data <- impute_knn_qc_neg_peaklist2

# Create vector containing block information
sample_names <- rownames(pca_data)
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

pca_data <- cbind(pca_data, block)
write.table(pca_data, file = "impute_knn_pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Now do PCA
pca_data <- read.table(file = "impute_knn_pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, shape = "block", colour = "block", main = 'PCA on negative QC data with imputed k-means for missing values', frame = TRUE, frame.type = 'norm')
ggsave("impute_knn_pca_data.png")


######################################
# Filter k-means imputed data by RSD #
######################################

# Prepare data for RSD analysis
rsd_qc_neg_peaklist <- impute_knn_qc_neg_peaklist2
rsd_qc_neg_peaklist <- t(rsd_qc_neg_peaklist)

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
  feature_vector_rsds <- apply(data, 1, function(data) {
    sd <- sd(data)
    mean <- mean(data)
    rsd <- (sd/mean) * 100
    return(rsd)
  })

  # Column bind feature_vector_means and feature_vector_standard_deviations
  data <- cbind(data, feature_vector_means, feature_vector_standard_deviations,feature_vector_rsds)
  # List peak IDs from rownames which are above 20% RSD and need removing from
  # feature vector matrix
  rows_for_deleting <- which(data[, "feature_vector_rsds"] > 20)
  data <- data[-rows_for_deleting, ]
  return(data)
}

rsd_qc_neg_peaklist <- filter_by_rsd(rsd_qc_neg_peaklist)

############################################################
# Do PCA to check effect of RSD filtering after imputation #
############################################################

pca_data <- rsd_qc_neg_peaklist

# Remove RSD calculation columns
rsd_cols <- colnames(pca_data) %in% c("feature_vector_means", "feature_vector_standard_deviations", "feature_vector_rsds")
pca_data <- pca_data[,!rsd_cols]
pca_data <- t(pca_data)

# Create vector containing block information
sample_names <- rownames(pca_data)
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

pca_data <- cbind(pca_data, block)
write.table(pca_data, file = "rsd_knn_neg_pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

pca_data <- read.table(file = "rsd_knn_neg_pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, shape = "block", colour = "block", main = 'PCA on RSD filtered, k-means imputed QC negative data', frame = TRUE, frame.type = 'norm')
ggsave("rsd_neg_qc_pca_data.png")


################################################
# Get feature vectors from samples and combine #
# with negative QCs                            #
################################################

pretreated_qc_sample_peaklist <- rsd_qc_neg_peaklist
pretreated_qc_sample_peaklist <- pretreated_qc_sample_peaklist[, !rsd_cols]
# Get peak idx of pretreated QC feature vectors
pretreated_idx <- rownames(pretreated_qc_sample_peaklist)

# Get QC and sample feature vectors from neg_peaklist
pretreated_sample_peaklist <- neg_peaklist[rownames(neg_peaklist) %in% pretreated_idx,]
# Get rid of QCs in pretreated_sample_peaklist
neg_sample_names <- meta_all[meta_all$type == "Sample", "file_name_neg"]
neg_sample_names <- as.character(neg_sample_names)
# Get sample peaklist using neg_sample_names
pretreated_sample_peaklist <- pretreated_sample_peaklist[,colnames(pretreated_sample_peaklist) %in% neg_sample_names]
# Combine QCs and samples from negative peaklist
pretreated_qc_sample_peaklist <- cbind(pretreated_qc_sample_peaklist, pretreated_sample_peaklist)


#########################################################################
# Do PCA on combined rsd_qc_neg_peaklist and pretreated_sample_peaklist #
#########################################################################

pca_data <- pretreated_qc_sample_peaklist
pca_data_qc_sample_names <- colnames(pretreated_qc_sample_peaklist)
pca_data_peak_idx <- rownames(pretreated_qc_sample_peaklist)
pca_data <- t(pca_data)

# Create vector containing block information
sample_names <- rownames(pca_data)
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

# Prepare QC and sample information for labelling data points
pca_meta_qc_sample <- integer(0)
qc_sample_names <- colnames(pretreated_qc_sample_peaklist)
for (i in 1:length(qc_sample_names)) {
  if (meta_all[which(meta_all[,"file_name_neg"]==qc_sample_names[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

pca_data <- cbind(pca_data, block, pca_meta_qc_sample)
# PCA cannot be performed on data with missing values
# Therefore remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Colour blocks in PCA plot
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:210]), data = pca_data, colour = "block", shape = "pca_meta_qc_sample", main = 'PCA on pretreated negative QCs with negative samples', frame = TRUE, frame.type = 'norm')
ggsave("rsd_neg_qc_get_samples_coloured_blocks.png")

# Colour QCs and samples in PCA plot
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:210]), data = pca_data, colour = "pca_meta_qc_sample", shape = "block", main = 'PCA on pretreated negative QCs with negative samples', frame = TRUE, frame.type = 'norm')
ggsave("rsd_neg_qc_get_samples.png")


############################################################
# Use k-means to impute missing values in negative samples #
############################################################

knn_sample_peaklist <- pretreated_sample_peaklist
# Transpose data
knn_sample_peaklist <- t(knn_sample_peaklist)

# Identify which column peak features have missing values
na_cols <- colnames(knn_sample_peaklist)[colSums(is.na(knn_sample_peaklist)) > 0]
# Perform K-means clustering using VIM
knn <- kNN(knn_sample_peaklist, variable = na_cols)

# Copy peak features from kNN result
knn_sample_peaklist2 <- knn[,1:506]
knn_sample_peaklist2 <- as.matrix(knn_sample_peaklist2)
# Copy rownames from impute_knn_data to rownames for knn
rownames(knn_sample_peaklist2) <- rownames(knn_sample_peaklist)
# Copy colnames from impute_knn_data to colnames for knn
colnames(knn_sample_peaklist2) <- colnames(knn_sample_peaklist)

# Combine QCs and knn samples from negative peaklist
knn_sample_peaklist2 <- t(knn_sample_peaklist2)
pretreated_qc_knn_sample_peaklist <- cbind(rsd_qc_neg_peaklist[,!rsd_cols], knn_sample_peaklist2)


#########################################################################
# Do PCA on combined rsd_qc_neg_peaklist and pretreated_sample_peaklist #
#########################################################################

pca_data <- pretreated_qc_knn_sample_peaklist
pca_data_qc_sample_names <- colnames(pretreated_qc_knn_sample_peaklist)
pca_data_peak_idx <- rownames(pretreated_qc_knn_sample_peaklist)
pca_data <- t(pca_data)

# Create vector containing block information
sample_names <- rownames(pca_data)
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

# Prepare QC and sample information for labelling data points
pca_meta_qc_sample <- integer(0)
qc_sample_names <- rownames(pca_data)
for (i in 1:length(qc_sample_names)) {
  if (meta_all[which(meta_all[,"file_name_neg"]==qc_sample_names[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

pca_data <- cbind(pca_data, block, pca_meta_qc_sample)
write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:506]), data = pca_data, colour = "pca_meta_qc_sample", shape = "block", main = 'PCA on pretreated negative QCs with knn negative samples', frame = TRUE, frame.type = 'norm')
ggsave("pretreated_qc_knn_neg_samples_pca.png")

pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:506]), data = pca_data, colour = "pca_meta_qc_sample", shape = "block", main = 'PCA on pretreated negative QCs with knn negative samples', frame = TRUE, frame.type = 'norm')
ggsave("pretreated_qc_knn_neg_samples_pca.png")


############################################################
# Use metaX data normalisation to correct within batch and #
# batch-to-batch variation in negative QCs and samples     #
############################################################

# Need to consider batch differences affecting peaks in QCs due to the
# 4 analytical blocks are removed by XCMS during peak alignment.

# Data for metaX signal correction needs to look like this:
head(signal_correction_data[,1:4])
## name batch01_QC01 batch01_QC02 batch01_QC03
## 1  78.02055      14023.0      13071.0      15270.0
## 2 452.00345      22455.0      10737.0      27397.0
## 3 138.96337       6635.4       8062.3       6294.6
## 4  73.53838      26493.0      26141.0      25944.0
## 5 385.12885      57625.0      56964.0      59045.0
## 6 237.02815     105490.0      90166.0      92315.0

signal_corr_data <- pretreated_qc_knn_sample_peaklist

# Move rownames into first column
signal_corr_data <- cbind(rownames(signal_corr_data), signal_corr_data)
# Check
head(signal_corr_data[ , 1:4])
##         GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_16r
## 2  "2"  "2214161.06684043" "2213167.30024468" "2200922.27885106"
## 11 "11" "80404.2911244445" "90185.9247093023" "80860.0580355555"
## 12 "12" "1052155.88762273" "1013478.89673721" "968455.803136363"
## 46 "46" "95204.1498043478" "95296.5965466666" "96498.8371822222"
## 47 "47" "37436017.503675"  "36930624.0750681" "37326859.8432458"
## 51 "51" "47954.2056"       "43746.2516363636" "53033.2722418605"

# Add colname for first column
col_names <- c("name", colnames(signal_corr_data)[-1])
colnames(signal_corr_data) <- col_names

# Save as file and reload again
write.table(signal_corr_data, file = "signal_corr_data.tab", sep ="\t", row.names = TRUE, col.names = TRUE)
signal_corr_data <- read.table(file = "signal_corr_data.tab", sep="\t")

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
# Save as file and reload again
write.table(sampleListFile, file = "sampleListFile.tab", sep ="\t", row.names = TRUE, col.names = TRUE)
sampleListFile <- read.table(file = "sampleListFile.tab", sep="\t")

# Do normalisation
para <- new("metaXpara")
pfile <- paste(datadir, "gsk/raw/esi_neg/netcdf/signal_corr_data.tab", sep="")
sfile <- paste(datadir, "gsk/raw/esi_neg/netcdf/sampleListFile.tab", sep="")
rawPeaks(para) <- read.delim(pfile, check.names = FALSE)
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
# para <- missingValueImpute(para)
para <- doQCRLSC(para, cpu=12)
plotQCRLSC(para$metaXpara)

# Download data
non_norm_neg_data <- getPeaksTable(para, valueID="value")
write.table(non_norm_neg_data, file = "non_norm_neg_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
norm_neg_data <- getPeaksTable(para, valueID="valueNorm")
write.table(norm_neg_data, file = "norm_neg_data.csv", sep =",", row.names = TRUE, col.names = TRUE)


#####################################
# Do PCA to check signal correction #
#####################################

head(rsd_knn_qc_neg[, 1:4])
##    GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_16r GSK_neg_block1_21r
## 2          2214161.07         2213167.30         2200922.28         2229604.08
## 11           80404.29           90185.92           80860.06           70636.49
## 12         1052155.89         1013478.90          968455.80          955786.38
## 46           95204.15           95296.60           96498.84           76137.69
## 47        37436017.50        36930624.08        37326859.84        34001468.40
## 51           47954.21           53991.18           53033.27           34404.22

head(peak_data)
##   ID             sample       value batch class order
## 1  2 GSK_neg_block1_09r  2214161.07     1    NA     1
## 2 11 GSK_neg_block1_09r    80404.29     1    NA     1
## 3 12 GSK_neg_block1_09r  1052155.89     1    NA     1
## 4 46 GSK_neg_block1_09r    95204.15     1    NA     1
## 5 47 GSK_neg_block1_09r 37436017.50     1    NA     1
## 6 51 GSK_neg_block1_09r    47954.21     1    NA     1

# Get signal corrected peak data
para <- transformation(para$metaXpara, valueID = "valueNorm")
norm_peak_data <- getPeaksTable(para, valueID = "valueNorm")
# Use sample column as rownames
rownames(norm_peak_data) <- norm_peak_data$sample
# Remove sample, class, batch and order columns
drop_cols <- c("sample", "class", "batch", "order")
norm_peak_data <- norm_peak_data[ , !(names(norm_peak_data) %in% drop_cols)]
# Transpose data
norm_peak_data <- t(norm_peak_data)
# Order data based on rownames (which are peak IDs)
norm_peak_data <- norm_peak_data[naturalorder(rownames(norm_peak_data)),]

pca_data <- norm_peak_data

# Prepare block information for labelling data points
block <- integer(0)
sample_names <- colnames(pca_data)
for (i in 1:length(qc_sample_names)) {
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

# Prepare QC and sample information for labelling data points
pca_meta_qc_sample <- integer(0)
for (i in 1:length(sample_names)) {
  if (meta_all[which(meta_all[,"file_name_neg"]==sample_names[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

regimens <- lapply(sample_names, function(sample_name) {
  regimen <- meta_all[which(meta_all[,"file_name_neg"] == sample_name), "Regimen"]
  if(is.na(regimen)) {
    regimen <- "QC"
  }
  regimen <- as.character(regimen)
})
regimens <- unlist(regimens)

# Transpose data
pca_data <- t(pca_data)
# Add block information, QC and sample names to PCA data
pca_data <- cbind(pca_data, block, pca_meta_qc_sample, regimens)

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")

# Highlight sample type and batch in PCA plot
autoplot(prcomp(pca_data[, 1:506]), data = pca_data, shape= "block", colour = 'pca_meta_qc_sample', main = 'PCA on normalised negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("norm_neg_qc_sample_pca_coloured_qc_sample.png")

# Highlight regimens in PCA plot
autoplot(prcomp(pca_data[, 1:506]), data = pca_data, colour = 'regimens', main = 'PCA on normalised negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("norm_neg_qc_sample_pca_coloured_regimens.png")

# Highlight blocks in PCA plot
autoplot(prcomp(pca_data[, 1:506]), data = pca_data, shape = "pca_meta_qc_sample", colour = 'block', main = 'PCA on normalised negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("norm_neg_qc_sample_pca_coloured_blocks.png")


#############################################################################
# Perform Canonical Variate Analysis on food and exercise regimens with QCs #
#############################################################################

ht(iris)
## Sepal.Length Sepal.Width Petal.Length Petal.Width   Species
## 1            5.1         3.5          1.4         0.2    setosa
## 2            4.9         3.0          1.4         0.2    setosa
## 3            4.7         3.2          1.3         0.2    setosa
## 4            4.6         3.1          1.5         0.2    setosa
## 5            5.0         3.6          1.4         0.2    setosa
## 6            5.4         3.9          1.7         0.4    setosa
## 7            4.6         3.4          1.4         0.3    setosa
## 8            5.0         3.4          1.5         0.2    setosa
## 9            4.4         2.9          1.4         0.2    setosa
## 10           4.9         3.1          1.5         0.1    setosa
## 141          6.7         3.1          5.6         2.4 virginica
## 142          6.9         3.1          5.1         2.3 virginica
## 143          5.8         2.7          5.1         1.9 virginica
## 144          6.8         3.2          5.9         2.3 virginica
## 145          6.7         3.3          5.7         2.5 virginica
## 146          6.7         3.0          5.2         2.3 virginica
## 147          6.3         2.5          5.0         1.9 virginica
## 148          6.5         3.0          5.2         2.0 virginica
## 149          6.2         3.4          5.4         2.3 virginica
## 150          5.9         3.0          5.1         1.8 virginica

vari <- iris[,1:4]
h(vari)
## Sepal.Length Sepal.Width Petal.Length Petal.Width
## 1          5.1         3.5          1.4         0.2
## 2          4.9         3.0          1.4         0.2
## 3          4.7         3.2          1.3         0.2
## 4          4.6         3.1          1.5         0.2
## 5          5.0         3.6          1.4         0.2
## 6          5.4         3.9          1.7         0.4

facto <- iris[,5]
h(facto)
## [1] setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica

pca_results <- prcomp(pca_data[, 1:506])

# Do CVA
groups <- as.factor(regimens)
cvall <- CVA(pca_results$x[, 1:5], groups)

# Plot CVA results
pdf('regimens-cva.pdf')
#png(file = "myplot.png", bg = "white")
plot(cvall$CVscores, col=groups, pch=as.numeric(groups), typ="n", asp=1,
xlab=paste("1st canonical axis", paste(round(cvall$Var[1, 2], 1), "%")),
ylab=paste("2nd canonical axis", paste(round(cvall$Var[2, 2], 1), "%")))
text(cvall$CVscores, as.character(groups), col=as.numeric(groups), cex=.7)

# add chull (merge groups)
for(jj in 1:length(levels(groups))){
  ii=levels(groups)[jj]
  kk=chull(cvall$CVscores[groups==ii,1:2])
  lines(cvall$CVscores[groups==ii,1][c(kk, kk[1])],
  cvall$CVscores[groups==ii, 2][c(kk, kk[1])], col=jj)
}

# add 80% ellipses
if (require(car)) {
for(ii in 1:length(levels(groups))){
  dataEllipse(cvall$CVscores[groups==levels(groups)[ii],1],
  cvall$CVscores[groups==levels(groups)[ii],2],
  add=TRUE, levels=.80, col=c(1:7)[ii])}
}

dev.off()


#########
# ROPLS #
#########

norm_peaks <- getPeaksTable(para, valueID="valueNorm")
# Create dataMatrix object
dataMatrix <- norm_peaks[c(1, 5:ncol(norm_peaks))]
rownames(dataMatrix) <- dataMatrix[, 1]
dataMatrix <- dataMatrix[ -c(1) ]

# Create sampleMetadata object which needs to look like this:
## > head(sampleMetadata)
## age   bmi gender
## HU_011  29 19.75      M
## HU_014  59 22.64      F
## HU_015  42 22.72      M
## HU_017  41 23.03      M
## HU_018  34 20.96      M
## HU_019  35 23.41      M

sampleMetadata <- norm_peaks[c(2:4)]
rownames(sampleMetadata) <- norm_peaks[, 1]

# variableMetadata object contains the chemical identity of the peaks
# Have not created this variableMetadata object

# Run PCA
ropls.norm.data.pca <- opls(dataMatrix)

# Save PCA plot
pdf('gsk_pca.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("overview", "outlier", "x-score", "x-loading"))
plot(ropls.norm.data.pca, typeVc = typeC, parDevNewL = FALSE)
dev.off()

# Check if there is an partitioning based on regimen class
regimenFc <- sampleMetadata[, "class"]
pdf('gsk_pca_regimen.pdf')
plot(ropls.norm.data.pca, typeVc = "x-score", parAsColFcVn = regimenFc, parEllipsesL = TRUE, parDevNewL = FALSE)
dev.off()

# Perform PLS and PLS-DA
gsk.plsda <- opls(dataMatrix, regimenFc, permI = 100, plotL = FALSE)
pdf('gsk_plsda_regimen.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("permutation", "overview", "outlier", "x-score"))
plot(gsk.plsda, typeVc = typeC, parDevNewL = FALSE)
dev.off()


#######
# PCA #
#######

para <- transformation(para$metaXpara, valueID = "valueNorm")
metaX::plotPCA(para, valueID = "valueNorm", scale = "pareto", center = TRUE, rmQC = FALSE)

##########
# PLS-DA #
##########

para <- transformation(para$metaXpara, valueID = "valueNorm")
plsdaPara <- new("plsDAPara")
plsdaPara@nperm <- 100
plsda.res <- runPLSDA(para = para, plsdaPara = plsdaPara, sample = c("NA", "A", "B", "C", "D"), valueID = "valueNorm")


#########################
# Normalise data to QCs #
#########################

# Get QC data
qc_names <- meta_all[meta_all[, "type"] == "QC", "file_name_neg"]
qc_names <- as.character(qc_names)
qc_data <- sig_corr_qc_sample_neg[ ,qc_names]

# Calculate average value for each QC peak
qc_means <- rowMeans(qc_data)

# Divide all data with QC data
qcnorm_qc_sample_neg <- sig_corr_qc_sample_neg[ , qc_sample_names]
# new_pca_data <- matrix(nrow = nrow(pca_data), ncol = ncol(pca_data))
# colnames(new_pca_data) <- colnames(pca_data)
# rownames(new_pca_data) <- rownames(pca_data)

# Test
qcnorm_qc_sample_neg <- cbind(qcnorm_qc_sample_neg, qc_means)
qcnorm_qc_sample_neg <- as.matrix(qcnorm_qc_sample_neg)
counter <- 0
ans <- apply(qcnorm_qc_sample_neg, 1, function(x) {
  counter <<- counter + 1
  # print(x)
  # print(paste("Length of x:", length(x)))
  qc_average <- as.numeric(x[length(x)])
  # print("##### End of matrix row #####")
  # for(i in 1:length(x)-1) {
  # print(x[length(x)])
  #   print(x[i] / x[length(x)])
  # }
  stuff <- lapply(x, function (y) {
    # print(paste("y:", y))
    # qc_average <- as.numeric(y[length(y)])
    # print(paste("qc_average: ", qc_average))
    new_value <- y / qc_average
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
  qcnorm_qc_sample_neg[counter,] <<- stuff
})

# Clean up
qcnorm_qc_sample_neg <- qcnorm_qc_sample_neg[, -ncol(qcnorm_qc_sample_neg)]

# Do PCA
pca_data <- qcnorm_qc_sample_neg

# Transpose data
pca_data <- t(pca_data)
# Add regimen information to PCA data
pca_data <- cbind(pca_data, regimen)

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = 'regimen', main = 'PCA on QC-normalised negative QC and sample data')
ggsave("pca_qc_normalised.png")


##################################
# Normalise data to time point 0 #
##################################

# Get all TP 0 data
tp0_sample_names <- meta_all[meta_all[, "Timepoint..PIMS."] == 0, "file_name_neg"]
tp0_sample_names <- as.character(tp0_sample_names)
tp0_sample_names <- na.omit(tp0_sample_names)
tp0_data <- qcnorm_qc_sample_neg[ ,tp0_sample_names]

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