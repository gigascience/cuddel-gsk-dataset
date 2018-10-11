# Title     : formatData.R
# Objective : Format data to create dataMatrix, sampleMetadata and
#             variableMetadata files
# Created by: peterli
# Created on: 11/10/2018

source("functions.R")

# Location of GSK data set
datadir = "/home/peter/"
# File path for positive files
pos_dir = paste0(datadir, "gsk/raw/esi_pos/netcdf")
# Output directory
output_path <- paste0(pos_dir, "/output")

# Create dataMatrix for positive-mode data
# Read in pre-treated positive mode data
pos_tp0_peaklist <- read.csv(paste(output_path, "pos_tp0_norm_plasma_peaklist.csv", sep="/"))
colnames <- colnames(pos_tp0_peaklist)
pos_cols <- grep('_pos_', colnames, value=TRUE)
# Sort out feature variables ID
dataMatrix <- paste0("P", rownames(pos_tp0_peaklist))
# Create dataMatrix for positive-mode data
pos_dataMatrix <- cbind(dataMatrix, pos_tp0_peaklist[, pos_cols])
# Sort out rownames
rownames(pos_dataMatrix) <- seq(1, nrow(pos_dataMatrix))

# Create variableMetadata for positive-mode data
# Define columns not to extract
meta_cols <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")
# Sort out feature variables ID
variableMetadata <- paste0("P", rownames(pos_tp0_peaklist))
# Read in pre-treated positive mode data
pos_xcms_peaklist <- read.csv(paste(output_path, "xcms_pos_peaklist.csv", sep="/"))
# Create variableMetadata matrix for positive-mode data
pos_variableMetadata <- cbind(variableMetadata, pos_xcms_peaklist[rownames(pos_tp0_peaklist), meta_cols])
# Sort out rownames
rownames(pos_variableMetadata) <- seq(1, nrow(pos_dataMatrix))

# Create sampleMetadata for positive-mode data
# Read in metadata from file
meta <- read.csv(paste0(datadir, "gsk/meta/meta_20180831.csv"))
meta <- meta[which(meta[, "file_name_pos"] %in% colnames), ]
meta_colnames <- colnames(meta)
req_cols <- c("file_name_pos", "sample_ID", "order", "type", "block", "Regimen", "food", "exercise")
meta <- meta[, req_cols]
# Correct column rownames
colnames(meta)[3] <- "injectionOrder"
colnames(meta)[4] <- "sampleType"
colnames(meta)[5] <- "batch"
colnames(meta)[6] <- "regimen"
mode <- replicate(nrow(meta), "pos")
meta <- cbind(meta, mode)
colnames(meta)[1] <- "sampleMetadata"
rownames(meta) <- seq(1, nrow(meta))
