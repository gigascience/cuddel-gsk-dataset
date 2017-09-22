library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(mice)
library(ggfortify)
library(metaX)
library(VIM)

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

#################################################
# Do PCA on unprocessed negative QC and samples #
#################################################

# Prepare neg_peaklist data
pca_data <- neg_peaklist[, file_name_neg]

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
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape= "block", colour = 'pca_meta_qc_sample', main = 'PCA on unprocessed negative QC and sample data')
ggsave("unprocessed_neg_qc_sample_pca.png")


###################################################
# Extract negative QC data from negative peaklist #
###################################################

# Get QC sample names
meta_all_qc_rows <- meta_all[, "type"] == "QC"
neg_qc_names <- meta_all[meta_all_qc_rows, "file_name_neg"]
neg_qc_names <- as.character(neg_qc_names)
# Subset QC data from negative peak list
qc_neg_peaklist <- neg_peaklist[,neg_qc_names]


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
autoplot(prcomp(pca_data[,1:409]), data = pca_data, colour= "block", main = 'PCA on unprocessed negative QCs')
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
del40_qc_neg_peaklist <- percent_na_qc_neg_peaklist[over_40_percent_rows,]


###############################################################
# Do PCA after removing features with 40% missing data values #
###############################################################

# Not required columns for PCA analysis
nr <- c("percent_nas", "idx")
nr <- colnames(del40_qc_neg_peaklist) %in% c("percent_nas", "idx")
# Prepare neg_peaklist data
pca_data <- del40_qc_neg_peaklist[,!nr]

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
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour= "block", main = 'PCA on negative QCs after peaks missing 40% values removed')
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
##################################################################################

# Replace NAs with k-means clustering
impute_knn_qc_neg_peaklist <- del40_qc_neg_peaklist
# Clean up columns
nr <- c("percent_nas", "idx")
nr <- colnames(del40_qc_neg_peaklist) %in% c("percent_nas", "idx")
impute_knn_qc_neg_peaklist <- impute_knn_qc_neg_peaklist[,!nr]
# Transpose data
impute_knn_qc_neg_peaklist <- t(impute_knn_qc_neg_peaklist)
# Identify which column peak features have missing values
na_cols <- colnames(impute_knn_qc_neg_peaklist)[colSums(is.na(impute_knn_qc_neg_peaklist)) > 0]
# Perform K-means clustering using VIM
knn <- kNN(impute_knn_qc_neg_peaklist, variable = na_cols)
# Copy peak features from kNN result
impute_knn_qc_neg_peaklist2 <- knn[,1:1704]
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
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = "block", main = 'PCA on negative QC data with imputed k-means for missing values')
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
  rows_for_deleting <- which(data[,"feature_vector_rsds"] > 20)
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
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = "block", main = 'PCA on RSD filtered, k-means imputed QC and sample negative data')
ggsave("rsd_neg_qc_pca_data.png")


#####################################################
# Get feature vectors from samples using QC vectors #
#####################################################

pretreated_qc_sample_peaklist <- rsd_qc_neg_peaklist
pretreated_qc_sample_peaklist <- pretreated_qc_sample_peaklist[,!rsd_cols]
# Get peak idx of pretreated QC feature vectors
pretreated_idx <- rownames(pretreated_qc_sample_peaklist)

# Get QC and sample feature vectors from neg_peaklist
pretreated_sample_peaklist <- neg_peaklist[rownames(neg_peaklist) %in% pretreated_idx,]
# Get rid of QCs in pretreated_sample_peaklist
neg_sample_names <- meta_all[meta_all$type == "Sample", "file_name_neg"]
neg_sample_names <- as.character(neg_sample_names)
# Get sample peaklist using neg_sample_names
pretreated_sample_peaklist <- pretreated_sample_peaklist[,colnames(pretreated_sample_peaklist) %in% neg_sample_names]


# Do PCA on combined rsd_qc_neg_peaklist and pretreated_sample_peaklist


############################################################
# Use metaX data normalisation to correct within batch and #
# batch-to-batch variation                                 #
############################################################

# Need to consider batch differences affecting peaks in QCs due to the
# 4 analytical blocks are removed by XCMS during peak alignment.

# Data for metaX signal correction needs to look like this:
> head(signal_correction_data[,1:4])
name batch01_QC01 batch01_QC02 batch01_QC03
1  78.02055      14023.0      13071.0      15270.0
2 452.00345      22455.0      10737.0      27397.0
3 138.96337       6635.4       8062.3       6294.6
4  73.53838      26493.0      26141.0      25944.0
5 385.12885      57625.0      56964.0      59045.0
6 237.02815     105490.0      90166.0      92315.0

# Remove unrequired columns
signal_corr_data <- subset(neg_peaklist, select = -c(1:12, 384:386))

# Remove peak rows if it contains missing values
# signal_corr_data <- signal_corr_data[ , colSums(is.na(signal_corr_data)) == 0]
signal_corr_data <- na.omit(signal_corr_data)

# Move rownames into first column
signal_corr_data <- cbind(rownames(signal_corr_data), signal_corr_data)
# Check
head(signal_corr_data[,1:4])
# Add colname for first column
col_names <- c("name", colnames(signal_corr_data)[-1])
colnames(signal_corr_data) <- col_names
# Reformat row names
rownames(signal_corr_data) <- rep(1:nrow(signal_corr_data))
# Save as file and reload again
write.table(signal_corr_data, file = "signal_corr_data.tab", sep ="\t", row.names = TRUE, col.names = TRUE)
signal_corr_data <- read.table(file = "signal_corr_data.tab", sep="\t")

# Create sample list file
## sample batch class order
## 1 batch01_QC01 1  <NA>     1
## 2 batch01_QC02 1  <NA>     2
sample <- colnames(signal_corr_data[-1])
batch <- subset(meta_all, file_name_neg %in% sample, "block")
class <- subset(meta_all, file_name_neg %in% sample, "Regimen")
order <- meta_all$order

sampleListFile <- cbind(sample, batch, class, meta_all$order)
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
res <- doQCRLSC(para, cpu=1)
plotQCRLSC(res$metaXpara)


#####################################
# Do PCA to check signal correction #
#####################################

> head(rsd_knn_qc_neg[,1:4])
GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_16r GSK_neg_block1_21r
2          2214161.07         2213167.30         2200922.28         2229604.08
11           80404.29           90185.92           80860.06           70636.49
12         1052155.89         1013478.90          968455.80          955786.38
46           95204.15           95296.60           96498.84           76137.69
47        37436017.50        36930624.08        37326859.84        34001468.40
51           47954.21           53991.18           53033.27           34404.22

> head(peak_data)
ID             sample       value batch class order
1  2 GSK_neg_block1_09r  2214161.07     1    NA     1
2 11 GSK_neg_block1_09r    80404.29     1    NA     1
3 12 GSK_neg_block1_09r  1052155.89     1    NA     1
4 46 GSK_neg_block1_09r    95204.15     1    NA     1
5 47 GSK_neg_block1_09r 37436017.50     1    NA     1
6 51 GSK_neg_block1_09r    47954.21     1    NA     1

# Get signal corrected peak data
peak_data <- res$metaXpara
peak_data <- peak_data@peaksData
peak_data <- peak_data[order(peak_data$sample),]
# Re-format data so that QCs and samples are lined up column by column
# if rowname = sample name 1 then copy to column 1 at the end
sample_names <- peak_data[, "sample"]
sample_names <- unique(sample_names)
sample_names <- as.character(sample_names)
test <- subset(peak_data, sample == "GSK_neg_block1_09r")
test <- peak_data[ which(peak_data$sample=='GSK_neg_block1_09r'), "value"]

# Create empty matrix
pca_data <- matrix(, nrow = 409, ncol = length(sample_names))

for (i in 1:length(sample_names)) {
  col_data <- peak_data[ which(peak_data$sample==sample_names[i]), "value"]
  pca_data[, i] <- col_data
}
colnames(pca_data) <- sample_names

# Create block metadata information for PCA display
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

# Create vector containing QC and sample information
pca_meta_qc <- integer(0)
for (i in 1:length(sample_names)) {
  if (meta_all[which(meta_all[,"file_name_neg"]==sample_names[i]), "type"] == 'QC') {
    pca_meta_qc[i] <- "QC"
  }
  else {
    pca_meta_qc[i] <- "Sample"
  }
}

# Transpose data
pca_data <- t(pca_data)
# Add block information to PCA data
pca_data <- cbind(pca_data, block, pca_meta_qc)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
#pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape= "block", colour = 'pca_meta_qc', main = 'PCA on signal corrected negative QC data')
ggsave("signal_corr_neg_qc_pca.png")

