# Title     : functions.R
# Objective : Provide functions for querying metadata about the GSK longitudinal
#             dataset
# Created by: peterli
# Created on: 14/8/2018

library(naturalsort)

# Location of GSK data set
datadir = "/home/peter/"
# File path for negative and positive files
neg_dir = paste(datadir, "gsk/raw/esi_neg/netcdf", sep="")
pos_dir = paste(datadir, "gsk/raw/esi_pos/netcdf", sep="")

######################################
# Read in metadata for data analysis #
######################################

# Read in and sort metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta_20180831.csv", sep=""))
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all_20180831.csv", sep=""))

#' getFileNameColumn
#'
#' Returns the column name containing the filenames for samples measured in
#' negative or positive ion mode.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return file_path. A column name for positive or negative sample names
getFileNameColumn <- function(mode="negative") {
    colname <- "file_name_neg"
    if (mode == "positive") {
        colname <- "file_name_pos"
    }
    return(colname)
}

#' getGSKFilePath
#'
#' Returns the file path for a given sample measured in negative or positive
#' mode.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return file_path. A vector containing block metadata.
getGSKFilePath <- function(gsk_filename) {
    file_path <- ""
    if (grepl("neg", gsk_filename)) {
        if (grepl("block1", gsk_filename)) {
            file_path <- paste(neg_dir, "/block1neg/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block2", gsk_filename)) {
            file_path <- paste(neg_dir, "/block2neg/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block3", gsk_filename)) {
            file_path <- paste(neg_dir, "/block3neg/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block4", gsk_filename)) {
            file_path <- paste(neg_dir, "/block4neg/", gsk_filename, ".cdf", sep = "")
        }
    }
    else if (grepl("pos", gsk_filename)) {
        if (grepl("block1", gsk_filename)) {
            file_path <- paste(pos_dir, "/block1pos/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block2", gsk_filename)) {
            file_path <- paste(pos_dir, "/block2pos/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block3", gsk_filename)) {
            file_path <- paste(pos_dir, "/block3pos/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block4", gsk_filename)) {
            file_path <- paste(pos_dir, "/block4pos/", gsk_filename, ".cdf", sep = "")
        }
    }
    return(file_path)
}

#' getAllGSKFilePaths
#'
#' Returns a vector containing the file paths for all negative or positive QC
#' and plasma samples.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return file_paths. A vector containing file paths for negative or positive
#' samples.
getAllGSKFilePaths <- function(mode="negative") {
    file_paths <- character()
    colname <- getFileNameColumn(mode)
    files <- meta[, colname]
    files <- as.character(files)
    for (i in 1:length(files)) {
        file_paths[i] <- getGSKFilePath(files[i])
    }
    return(file_paths)
}

#' getBlockMetadata
#'
#' Returns the block information assoicated with a vector of sample names.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return block_metadata. A vector containing block metadata.
getBlockMetadata <- function(sample_names, mode="negative") {
    block_meta <- character(0)

    colname <- getFileNameColumn(mode)
    for (i in 1:length(sample_names)) {
        block_meta[i] <- meta[which(meta[, colname]==sample_names[i]), "block"]
    }
    return(block_meta)
}

#' getQCSampleMetadata
#'
#' Returns the sample type information assoicated with a vector of sample names.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return sample_meta. A vector containing sample type metadata.
getQCSampleMetadata <- function(samples, mode="negative") {
    colname <- getFileNameColumn(mode)
    sample_meta <- character(0)
    for (i in 1:length(samples)) {
        if (meta[which(meta[, colname]==samples[i]), "type"] == 'QC') {
            type <- meta[which(meta[, colname]==samples[i]), "type"]
            sample_meta[i] <- "QC"
        }
        else {
            sample_meta[i] <- "Sample"
        }
    }
    return(sample_meta)
}

#' getRegimenMetadata
#'
#' Returns the regimens assoicated with a vector of sample names.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return sample_meta. A vector containing regimen information.
getRegimenMetadata <- function(samples, mode="negative") {
    col_name <- getFileNameColumn(mode)

    sample_meta <- character(0)
    for (i in 1:length(samples)) {
        regimen <- meta[which(meta[, col_name]==samples[i]), "Regimen"]
        regimen <- as.character(regimen)
        if (regimen=="X") {
            sample_meta[i] <- "QC"
        }
        else {
            sample_meta[i] <- regimen
        }
    }
    return(sample_meta)
}

#' getQCSampleNames
#'
#' Returns all QC sample names which were analysed in negative or positive
#' ionization mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return qc_names. A vector containing negative or positive QC sample names.
getQCSampleNames <- function(mode="negative") {
    col_name <- getFileNameColumn(mode)

    qc_names <- meta[meta[, "type"] == 'QC', col_name]
    qc_names <- as.character(qc_names)
    return(qc_names)
}

#' getNumberOfQCs
#'
#' Returns the total number of QCs measured in negative or positive ionization
#' mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer representing the number negative or positive QCs.
getNumberOfQCs <- function(mode="negative") {
    qc_names <- getQCSampleNames(mode)
    return(length(qc_names))
}

#' getPlasmaSampleNames
#'
#' Returns a vector of plasma sample names which have been measured in negative
#' or positive ionization mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return plasma_names. A vector containing a list of negative or positive
#' plasma sample names.
getPlasmaSampleNames <- function(mode="negative") {
    col_name <- getFileNameColumn(mode)
    plasma_names <- meta[meta[, "type"] == 'Sample', col_name]
    plasma_names <- as.character(plasma_names)
    return(plasma_names)
}

#' getPlasmaSampleNamesByBlock
#'
#' Returns a vector of plasma sample names for a given block which have been
#' measured in negative or positive ionization mode in the GSK longitudinal data
#' set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return plasma_names. A vector containing a list of negative or positive
#' plasma sample names for a given block.
getPlasmaSampleNamesByBlock <- function(mode="negative", block="1") {
    col_name <- getFileNameColumn(mode)
    # Get data associated for specified block
    block_data <- meta[meta[, "block"] == block, ]
    # Get plasma sample names from block_data
    plasma_names <- block_data[block_data[, "type"] == 'Sample', col_name]
    plasma_names <- as.character(plasma_names)
    return(plasma_names)
}

#' getAllSampleNamesByBlock
#'
#' Returns a vector of all QC and plasma sample names for a given block which
#' have been measured in negative or positive ionization mode in the GSK
#' longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return plasma_names. A vector containing a list of all QC and plasma sample
#' names for a given block.
getAllSampleNamesByBlock <- function(mode="negative", block="1") {
    col_name <- getFileNameColumn(mode)
    # Get all sample names from block_data
    all_sample_names <- meta[meta[, "block"] == block, col_name]
    all_sample_names <- as.character(all_sample_names)
    return(all_sample_names)
}

#' getNumberOfPlasmaSamples
#'
#' Returns an integer showing the number of plasma samples which have been
#' measured in negative or positive ionization mode in the GSK longitudinal data
#' set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer showing number of negative or positive plasma sample names.
getNumberOfPlasmaSamples <- function(mode="negative") {
    plasma_names <- getPlasmaSampleNames(mode)
    return(length(plasma_names))
}

#' getTP0PlasmaSampleNames
#'
#' Returns an integer showing the number of plasma sample names which have been
#' measured at Time Point 0 in negative or positive ionization mode in the GSK
#' longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return tp0_sample_names. An integer showing number of positive or negative
#' plasma sample names.
getTP0PlasmaSampleNames <- function(mode="negative") {
    col_name <- getFileNameColumn(mode)
    tp0_sample_names <- meta[meta[, "Timepoint"] == 0, col_name]
    tp0_sample_names <- as.character(tp0_sample_names)
    tp0_sample_names <- na.omit(tp0_sample_names)
    return(tp0_sample_names)
}
