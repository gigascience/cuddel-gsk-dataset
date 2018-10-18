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

# Read in data
pos_dataMatrix <- read.table(file=paste(output_path, "dataMatrix.tsv", sep="/"), header=TRUE, row.names=1)
rownames(pos_dataMatrix) <- pos_dataMatrix[,1]
pos_dataMatrix <- pos_dataMatrix[, 2:length(pos_dataMatrix)]
pos_variableMetadata <- read.table(file=paste(output_path, "variableMetadata.tsv", sep="/"), header=TRUE, row.names=1)
pos_sampleMetadata <- read.table(file=paste(output_path, "sampleMetadata.tsv", sep="/"), header=TRUE, row.names=1)

#' Calculate means for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of meanvalues for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
calcMean <- function(x, fac) {
    tapply(x, fac, mean)
}

#' Calculate means for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of mean values calculated for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
calcStdDev <- function(x, fac) {
    tapply(x, fac, sd)
}

#' Calculate t-test for multiple groups
#'
#' @param x A vector of values for mean calculation
#' @param tp5RegARegB_regimen_fac A vector of values representing factors associated with x vector
#' @return The matrix of mean values calculated for each factor represented in tp5RegARegB_regimen_fac
#' @examples
#' calcMean(10, 1)
performTTest <- function(x, fac) {
    # Need to split data based on the 2 factor groups
    grp_data <- split(x, fac)
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

