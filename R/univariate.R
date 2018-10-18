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

#' Calculate means for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of meanvalues for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
calcMean <- function(x, tp5RegARegB_regimen_fac) {
    tapply(x, tp5RegARegB_regimen_fac, mean)
}

#' Calculate means for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of mean values calculated for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
calcStdDev <- function(x, tp5RegARegB_regimen_fac) {
    tapply(x, tp5RegARegB_regimen_fac, sd)
}

#' Calculate t-test for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of mean values calculated for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
performTTest <- function(x, tp5RegARegB_regimen_fac) {
    # Need to split data based on the 2 factor groups
    grp_data <- split(x, tp5RegARegB_regimen_fac)
    ttest <- t.test(grp_data[[1]], grp_data[[2]])
    return(ttest$p.value)
}

# Get timepoints
tp_fac <- factor(pos_sampleMetadata[, "timepoint"])
tp_lev <- levels(tp_fac)

# Extract sample IDs data based on timepoint 5
tp5_sample_meta <- pos_sampleMetadata[pos_sampleMetadata[, "timepoint"] == "5", ]
# Extract sample IDs for regimen A and B data from tp5_sample_meta
tp5RegARegB_sample_meta <- tp5_sample_meta[tp5_sample_meta[, "regimen"] == "A" | tp5_sample_meta[, "regimen"] == "B", ]
tp5RegARegB_regimen_fac <- factor(tp5RegARegB_sample_meta[, "regimen"])

t_pos_dataMatrix <- t(pos_dataMatrix)

# Extract peak feature values for Reg A and B
data <- t_pos_dataMatrix[rownames(tp5RegARegB_sample_meta), ]

# Call function to do t-tests
if(nrow(data) != length(tp5RegARegB_regimen_fac)) {
    stop("The number of rows in data needs to equal length of fac!!")
} else {
    means <- apply(data, 2, calcMean, fac=tp5RegARegB_regimen_fac)
    stdevs <- apply(data, 2, calcStdDev, fac=tp5RegARegB_regimen_fac)
    pvalues <- apply(data, 2, performTTest, fac=tp5RegARegB_regimen_fac)
}

