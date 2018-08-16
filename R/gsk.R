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
head(signal_correction_data[, 1:4])
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