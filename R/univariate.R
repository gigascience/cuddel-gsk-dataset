# Title     : univariate.R
# Objective : Analyse GSK dataset with univariate statistics
# Created by: peterli
# Created on: 9/10/2018

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
pos_dataMatrix <- cbind(pos_tp0_peaklist[, pos_cols])
# Sort out rownames
rownames(pos_dataMatrix) <- dataMatrix

# Output data
write.table(
    pos_dataMatrix,
    file="dataMatrix.tsv",
    sep="\t",
    row.names=TRUE,
    col.names=TRUE)

# Create variableMetadata for positive-mode data
# Define columns not to extract
meta_cols <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")
# Sort out feature variables ID
variableMetadata <- paste0("P", rownames(pos_tp0_peaklist))
# Read in pre-treated positive mode data
pos_xcms_peaklist <- read.csv(paste(output_path, "xcms_pos_peaklist.csv", sep="/"))
# Create variableMetadata matrix for positive-mode data
pos_variableMetadata <- cbind(pos_xcms_peaklist[rownames(pos_tp0_peaklist), meta_cols])
# Sort out rownames
rownames(pos_variableMetadata) <- variableMetadata
# Output data
write.table(
pos_variableMetadata,
    file="variableMetadata.tsv",
    sep="\t",
    row.names=TRUE,
    col.names=TRUE)

# Create sampleMetadata for positive-mode data
# Read in metadata from file
meta <- read.csv(paste0(datadir, "gsk/meta/meta_20180831.csv"))
meta <- meta[which(meta[, "file_name_pos"] %in% colnames), ]
meta_colnames <- colnames(meta)
req_cols <- c("sample_ID", "order", "type", "block", "Regimen", "food", "exercise", "Timepoint")
pos_sampleMetadata <- meta[, req_cols]
# Correct column rownames
colnames(pos_sampleMetadata)[2] <- "injectionOrder"
colnames(pos_sampleMetadata)[3] <- "sampleType"
colnames(pos_sampleMetadata)[4] <- "batch"
colnames(pos_sampleMetadata)[5] <- "regimen"
colnames(pos_sampleMetadata)[8] <- "timepoint"
mode <- replicate(nrow(pos_sampleMetadata), "pos")
pos_sampleMetadata <- cbind(pos_sampleMetadata, mode)
colnames(meta)[1] <- "sampleMetadata"
rownames(pos_sampleMetadata) <- meta[,"file_name_pos"]
# Output data
write.table(
    pos_sampleMetadata,
    file="sampleMetadata.tsv",
    sep="\t",
    row.names=TRUE,
    col.names=TRUE)

# Compare regimens
facFcVn <- factor(pos_sampleMetadata[, "regimen"])
facLevVc <- levels(facFcVn)

# Get timepoint
TPfac <- factor(pos_sampleMetadata[, "timepoint"])
TPlev <- levels(TPfac)

# Use tapply to analyse data
t_pos_dataMatrix <- t(pos_dataMatrix)
batchfac <- factor(pos_sampleMetadata[,"batch"])
tapply(t_pos_dataMatrix[,1], batchfac, mean)

# Use tapply from within a function
doStuff <- function(x) {
    tapply(x, batchfac, mean)
}

# Use the above function to run through all the peak features
t_pos_dataMatrix <- t(pos_dataMatrix)
staVn <- apply(t_pos_dataMatrix, 2, doStuff)

# Extract sample IDs data based on timepoint 5
tp5dat <- pos_sampleMetadata[pos_sampleMetadata[, "timepoint"] == "5", ]
# Extract sample IDs for regimen A and B data from tp5dat
tp5RegARegBdat <- tp5dat[tp5dat[, "regimen"] == "A" | tp5dat[, "regimen"] == "B", ]
fac <- factor(tp5RegARegBdat[, "regimen"])

#' Calculate means for multiple groups
#'
#' @param x A vector
#' @param fac A factor
#' @return The matrix of means for the groups.
#' @examples
#' add(1, 1)
#' add(10, 1)
doStuff2 <- function(x, fac) {
    tapply(x, fac, mean)
}
# Extract peak feature values for Reg A and B
data <- t_pos_dataMatrix[rownames(tp5RegARegBdat), ]
# Check length of data
print(paste("length of data", length(data)))
print(paste("length of factor", length(fac)))
# Call function
new2 <- apply(data, 2, doStuff2, fac=fac)

