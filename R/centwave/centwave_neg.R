# Title     : centwave_neg.R
# Objective : TODO
# Created by: peterli
# Created on: 16/11/2018

library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)

# Helper functions for selecting data files
source('functions.R')

#### Data import ####

# Location of GSK data set
datadir = "/Users/peterli"
# File path for negative files
neg_dir = paste(datadir, "gsk/raw/esi_neg/netcdf", sep="")
# Output directory
output_path <- paste(neg_dir, "/output", sep="")

######################################
# Read in metadata for data analysis #
######################################

# Get list of negative file paths
neg_file_paths <- getAllGSKFilePaths(mode="negative")
# Get list of negative QC and plasma sample names
qc_plasma_samples <- getAllQCPlasmaSampleNames(mode="negative")
# Get list of negative QC sample names
qc_samples <- getQCSampleNames(mode="negative")
# Get file paths for QC sample names
qc_sample_file_paths <- lapply(qc_samples, getGSKFilePath)
qc_sample_file_paths <- unlist(qc_sample_file_paths)

test_qcs <- getQCSampleNamesByBlock(block="1")[1:2]
test_qcs <- append(test_qcs, getQCSampleNamesByBlock(block="2")[1:2])
test_qcs <- append(test_qcs, getQCSampleNamesByBlock(block="3")[1:2])
test_qcs <- append(test_qcs, getQCSampleNamesByBlock(block="4")[1:2])

test_qcs_paths <- lapply(test_qcs, getGSKFilePath)
test_qcs_paths <- unlist(test_qcs_paths)

# Create a phenodata data.frame to describe experimental setup
pd <- data.frame(sample_name = sub(basename(test_qcs_paths), pattern = ".cdf", replacement = "", fixed = TRUE),
    sample_group = c(rep("block1", 2), rep("block2", 2), rep("block3", 2), rep("block4", 2)),
    stringsAsFactors = FALSE)

# Load raw data as OnDiskMSnExp object
raw_data <- readMSData(files = test_qcs_paths, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

#### Initial data inspection ####

# Look at retention time values
head(rtime(raw_data))

# Organise mx values by file
mzs <- mz(raw_data)
mzs_by_file <- split(mzs, f = fromFile(raw_data))  # Split the list by file
length(mzs_by_file)

## Get the base peak chromatograms. This reads data from the files.
bpis <- chromatogram(raw_data, aggregationFun = "max")
## Define colors for the two groups
group_colors <- paste0(brewer.pal(4, "Set1")[1:4], "60")
names(group_colors) <- c("block1", "block2", "block3", "block4")

## Plot all chromatograms - this is the total ion chromatogram
pdf("tic_plot.pdf")
plot(bpis, col = group_colors[raw_data$sample_group])
dev.off()


######################################
# Apply XCMS onto QC and sample data #
######################################

