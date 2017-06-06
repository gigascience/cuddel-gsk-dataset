library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(specmine)
library(impute)

# Read metadata file
meta_all <- read.csv("/home/peter/gsk/meta/meta_all.csv")
# Extract QCs
qc_meta_all <- subset(meta_all, type == "QC")
# Sort QCs by block then order
qc_meta_all <- qc_meta_all [ with(qc_meta_all, order(block, order)), ]
# Batch differences which might affect peaks in QCs due to the 4
# analytical blocks are removed by XCMS during peak alignment.
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


# Get list of neg qc filenames
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

# Get peak list df for neg QCs with index
# neg_qc_peaklist <- data.frame(neg_peaklist[1:12], neg_peaklist[neg_qc_fnames])
# Get peak list df for neg QCs
neg_qc_peaklist <- neg_peaklist[neg_qc_fnames]

# There should be 84 QCs
length(neg_qc_fnames)
[1] 84

# To get only QC data for row 1
# The 13th column is the first QC column
# neg_qc_peaklist[1,13:ncol(neg_qc_peaklist)]
neg_qc_peaklist[1, 1:ncol(neg_qc_peaklist)]

# Get number of NA in row 1
# sum(is.na(neg_qc_peaklist[1,13:ncol(neg_qc_peaklist)]))
sum(is.na(neg_qc_peaklist[1, 1:ncol(neg_qc_peaklist)]))

# Calculate percentage missing values in row 4
# 100*(sum(is.na(neg_qc_peaklist[4, 1:ncol(neg_qc_peaklist)]))/ncol(neg_qc_peaklist[ ,1:ncol(neg_qc_peaklist)]))
100*(sum(is.na(neg_qc_peaklist[4, 1:ncol(neg_qc_peaklist)]))/ncol(neg_qc_peaklist))

# Calculate percentage nas for each peak (row)
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

#################################################
# Get rid of rows with 40% or more missing values
#################################################
remove_bad_feature_vectors <- function(peak_list, neg_qc_peaklist, percentage_threshold = 40){
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

pretreated_peak_list1 <- remove_bad_feature_vectors(neg_qc_peaklist, percentage_threshold=40)

# Now need to fill in missing peak values
# Need to impute value for NA row by row since each row represent a specific
# peak
# Use specmine
# Need to create a list object which has a data object that is the peaklist containing IDs as rownames and QC names are column names
gsk_qc_neg <- list(data = neg_qc_peaklist[, 14:ncol(neg_qc_peaklist)], type = "ms-spectra", description = "GSK pooled QC")

impute_mean_gsk_qc_neg <- impute_nas_mean(gsk_qc_neg)
