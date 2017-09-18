library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(mice)
# library(impute)
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


###############################################################
# Calculate percentage of missing values for each feature row #
###############################################################

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

percent_nas <- rep(0, nrow(neg_peaklist))
percent_nas <- calculate_percentage_na(neg_peaklist)

# Paste percent_nas into negative peak list data frame
neg_peaklist <- cbind(neg_peaklist[1:12], percent_nas, neg_peaklist[file_name_neg])

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
del40_neg_peaklist <- remove_bad_feature_vectors(neg_peaklist)
del40_neg_peaklist <- as.matrix(as.data.frame(lapply(del40_neg_peaklist, as.numeric)))

###############################################################
# Do PCA after removing features with 40% missing data values #
###############################################################

# Prepare neg_peaklist data
pca_data <- del40_neg_peaklist[, file_name_neg]

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

# Prepare QC sample information for labelling data points
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
# Add block information to PCA data
pca_data <- cbind(pca_data, block, pca_meta_qc_sample)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape= "block", colour = 'pca_meta_qc_sample', main = 'PCA on negative QC and sample data after peaks missing 40% values removed')
ggsave("40_percent_filtered_neg_qc_sample_pca.png")

####################################################
# Impute means to replace remaining missing values #
####################################################

impute_mean_data <- del40_neg_peaklist[,14:ncol(del40_neg_peaklist)]
na_rows <- which(is.na(impute_mean_data), arr.ind=TRUE)
impute_mean_data[na_rows] <- rowMeans(impute_mean_data, na.rm=TRUE)[na_rows[,1]]


################################################
# Do PCA plot to check effect of imputed means #
################################################
impute_mean_pca_data <- impute_mean_data
rownames(impute_mean_pca_data) <- del40_neg_peaklist[,"idx"]
impute_mean_pca_data <- t(impute_mean_pca_data)

# Create vector containing block information
sample_names <- rownames(impute_mean_pca_data)
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

# Prepare QC sample information for labelling data points
pca_meta_qc_sample <- integer(0)
for (i in 1:length(rownames(impute_mean_pca_data))) {
  if (meta_all[which(meta_all[,"file_name_neg"]==rownames(impute_mean_pca_data)[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

impute_mean_pca_data <- cbind(impute_mean_pca_data, block, pca_meta_qc_sample)
write.table(impute_mean_pca_data, file = "impute_mean_pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Now do PCA
impute_mean_pca_data <- read.table(file = "impute_mean_pca_data.csv", sep=",")
autoplot(prcomp(impute_mean_pca_data[,1:1637]), data = impute_mean_pca_data, shape= "block", colour = "pca_meta_qc_sample", main = 'PCA on negative QC data with imputed means for missing values')
ggsave("impute_mean_pca_data.png")

##################################################################################
# Impute k-nearest neighbours to replace remaining missing values using specmine #
##################################################################################

# Replace NAs with k-means clustering
impute_knn_data <- del40_neg_peaklist[,14:ncol(del40_neg_peaklist)]
# Transpose data
impute_knn_data <- t(impute_knn_data)
# Copy rownames from del40_neg_peaklist to impute_knn_data column names
colnames(impute_knn_data) <- del40_neg_peaklist[,"idx"]
# Identify which column peak features have missing values
na_cols <- colnames(impute_knn_data)[colSums(is.na(impute_knn_data)) > 0]
# Perform K-means clustering using VIM
knn <- kNN(impute_knn_data, variable = na_cols)

na_rows <- which(is.na(impute_mean_data), arr.ind=TRUE)
impute_knn_data[na_rows] <- rowMeans(impute_mean_data, na.rm=TRUE)[na_rows[,1]]


#####################################################
# Do PCA plot to check effect of imputed knn values #
#####################################################
impute_mean_pca_data <- impute_mean_data
rownames(impute_mean_pca_data) <- del40_neg_peaklist[,"idx"]
impute_mean_pca_data <- t(impute_mean_pca_data)

# Create vector containing block information
sample_names <- colnames(del40_neg_peaklist)
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

# Prepare QC sample information for labelling data points
pca_meta_qc_sample <- integer(0)
for (i in 1:length(rownames(impute_mean_pca_data))) {
  if (meta_all[which(meta_all[,"file_name_neg"]==rownames(impute_mean_pca_data)[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

impute_mean_pca_data <- cbind(impute_mean_pca_data, block, pca_meta_qc_sample)
write.table(impute_mean_pca_data, file = "impute_mean_pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Now do PCA
impute_mean_pca_data <- read.table(file = "impute_mean_pca_data.csv", sep=",")
autoplot(prcomp(impute_mean_pca_data[,1:1648]), data = impute_mean_pca_data, shape= "block", colour = "pca_meta_qc_sample", main = 'PCA on negative QC data with imputed means for missing values')
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

