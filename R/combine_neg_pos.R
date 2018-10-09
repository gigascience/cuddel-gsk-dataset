# Title     : combine_neg_pos.R
# Objective : Combine pre-treated negative and positive mode datasets
# Created by: peterli
# Created on: 09/10/2018

# Location of GSK data set
datadir <- "/home/peter/"
# File path for negative files
neg_dir <- paste0(datadir, "gsk/raw/esi_neg/netcdf")
# File path for positive files
pos_dir <- paste0(datadir, "gsk/raw/esi_pos/netcdf")
# Set working directory for combined data
combined_dir <- paste0(datadir, "gsk/combined")
setwd(combined_dir)
# Load data files
neg_data <- read.csv(paste0(neg_dir, "/output/neg_tp0_norm_plasma_peaklist.csv"))
pos_data <- read.csv(paste0(pos_dir, "/output/pos_tp0_norm_plasma_peaklist.csv"))

# Create peak identifiers for negative data
PeakId <- rownames(neg_data)
PeakId <- lapply(PeakId, function(x){
    paste0(x, "_neg")
})
PeakId <- unlist(PeakId)
neg_data <- cbind(PeakId, neg_data)

# Clean up sample names for negative data
colnames <- colnames(neg_data)
colnames <- lapply(colnames, function(x){
    # If column name contains _neg_ or r then remove it
    x <- gsub("neg_", "", x)
    x <-gsub("r", "", x)
})
colnames <- unlist(colnames)
colnames(neg_data) <- colnames

# Do the same for positive data
PeakId <- rownames(pos_data)
PeakId <- lapply(PeakId, function(x){
    paste0(x, "_pos")
})
PeakId <- unlist(PeakId)
pos_data <- cbind(PeakId, pos_data)

# Clean up sample names for positive data
colnames <- colnames(pos_data)
colnames <- lapply(colnames, function(x){
    # If column name contains _neg_ or r then remove it
    x <- gsub("pos_", "", x)
    x <-gsub("r", "", x)
})
colnames <- unlist(colnames)
colnames(pos_data) <- colnames

# Combine pos and neg data
combined <- rbind(neg_data, pos_data)
# Create unique identifers for features in combined dataset
rownames(combined) <- c(1:nrow(combined))

# Output data
write.table(combined,
    file="combined_neg_pos.csv",
    sep=",",
    row.names=TRUE,
    col.names=TRUE)
