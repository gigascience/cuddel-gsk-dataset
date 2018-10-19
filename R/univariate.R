# Title     : univariate.R
# Objective : Analyse GSK dataset with univariate statistics
# Created by: peterli
# Created on: 9/10/2018

library(ggfortify)

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

# Get timepoints as a factor
tp_fac <- factor(pos_sampleMetadata[, "timepoint"])
# Get the levels of the timepoint factor
tp_levs <- levels(tp_fac)
# Remove -1 timepoint factor level as this corresponds to QCs
tp_levs <- tp_levs[2:length(tp_levs)]

# Transpose
t_pos_dataMatrix <- t(pos_dataMatrix)

# Create matrices to hold TP data for regimens A and B - 23 rows * 8 TP levels
tp_RegA_means <- matrix(0, nrow=23, ncol=8)
colnames(tp_RegA_means) <- c("TP0", "TP1", "TP4", "TP5", "TP8", "TP12", "TP16", "TP24")
rownames(tp_RegA_means) <- colnames(t_pos_dataMatrix)

tp_RegA_stdevs <- matrix(0, 23, 8)
colnames(tp_RegA_stdevs) <- c("TP0", "TP1", "TP4", "TP5", "TP8", "TP12", "TP16", "TP24")
rownames(tp_RegA_stdevs) <- colnames(t_pos_dataMatrix)

tp_RegB_means <- matrix(0, 23, 8)
colnames(tp_RegB_means) <- c("TP0", "TP1", "TP4", "TP5", "TP8", "TP12", "TP16", "TP24")
rownames(tp_RegB_means) <- colnames(t_pos_dataMatrix)

tp_RegB_stdevs <- matrix(0, 23, 8)
colnames(tp_RegB_stdevs) <- c("TP0", "TP1", "TP4", "TP5", "TP8", "TP12", "TP16", "TP24")
rownames(tp_RegB_stdevs) <- colnames(t_pos_dataMatrix)

tp_pvalues <- matrix(0, 23, 8)
colnames(tp_pvalues) <- c("TP0", "TP1", "TP4", "TP5", "TP8", "TP12", "TP16", "TP24")
rownames(tp_pvalues) <- colnames(t_pos_dataMatrix)

# Loop over timepoint levels
for(i in 1:length(tp_levs)) {
    # Extract sample IDs data based on timepoint 5
    tp_sample_meta <- pos_sampleMetadata[pos_sampleMetadata[, "timepoint"] == tp_levs[i], ]
    # Get regimen factor
    # Extract sample IDs for regimen A and B data from tp_sample_meta
    tp_reg_sample_meta <- tp_sample_meta[tp_sample_meta[, "regimen"] == "A" | tp_sample_meta[, "regimen"] == "B", ]
    # Get regimen factor
    regimen_fac <- factor(tp_reg_sample_meta[, "regimen"])

    # Extract peak feature values for Reg A and B
    data <- t_pos_dataMatrix[rownames(tp_reg_sample_meta), ]

    # Call function to do t-tests
    if(nrow(data) != length(regimen_fac)) {
        stop("The number of rows in data needs to equal length of fac!!")
    } else {
        means <- apply(data, 2, calcMean, fac=regimen_fac)
        t_means <- t(means)
        tp_RegA_means[, i] <- t_means[, "A"]
        tp_RegB_means[, i] <- t_means[, "B"]

        stdevs <- apply(data, 2, calcStdDev, fac=regimen_fac)
        t_stdevs <- t(stdevs)
        tp_RegA_stdevs[, i] <- t_stdevs[, "A"]
        tp_RegB_stdevs[, i] <- t_stdevs[, "B"]

        pvalues <- apply(data, 2, performTTest, fac=regimen_fac)
        t_pvalues <- t(pvalues)
        tp_pvalues[, i] <- t_pvalues
    }
}

t_tp_RegA_means <- t(tp_RegA_means)

# Basic line plot with points
ggplot(data=t_tp_RegA_means, aes(x=rownames(t_tp_RegA_means), y=t_tp_RegA_means[,1], group=1)) +
    geom_line() +
    geom_point()
# Change the line type
ggplot(data=df, aes(x=P210, y=len, group=1)) +
    geom_line(linetype = "dashed")+
    geom_point()
# Change the color
ggplot(data=df, aes(x=dose, y=len, group=1)) +
    geom_line(color="red")+
    geom_point()

ggsave("mtcars.pdf")