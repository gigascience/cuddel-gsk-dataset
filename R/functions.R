# Title     : TODO
# Objective : TODO
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
meta <- read.csv(paste(datadir, "gsk/meta/meta.csv", sep=""))
meta_neg <- meta[naturalorder(meta$file_name_neg), ]
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all.csv", sep=""))
meta_all_neg <- meta_all[naturalorder(meta_all$file_name_neg), ]

#' getAllNegGSKFilePaths
#'
#' Returns a vector containing all the file paths for negative samples.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return neg_qc_names. A vector containing file paths for negative samples.
getAllNegGSKFilePaths <- function() {
    neg_files <- meta_all[, "file_name_neg"]
    neg_files <- as.character(neg_files)
    neg_file_paths <- character()
    for (i in 1:length(neg_files)) {
        if (grepl("block1", neg_files[i])) {
            neg_file_path <- paste(neg_dir, "/block1neg/", neg_files[i], ".cdf", sep = "")
            neg_file_paths[i] <- neg_file_path
        }
        else if (grepl("block2", neg_files[i])) {
            neg_file_path <- paste(neg_dir, "/block2neg/", neg_files[i], ".cdf", sep = "")
            neg_file_paths[i] <- neg_file_path
        }
        else if (grepl("block3", neg_files[i])) {
            neg_file_path <- paste(neg_dir, "/block3neg/", neg_files[i], ".cdf", sep = "")
            neg_file_paths[i] <- neg_file_path
        }
        else if (grepl("block4", neg_files[i])) {
            neg_file_path <- paste(neg_dir, "/block4neg/", neg_files[i], ".cdf", sep = "")
            neg_file_paths[i] <- neg_file_path
        }
        else {
            neg_file_paths[i] <- ""
        }
    }
    return(neg_file_paths)
}

#' getAllPosGSKFilePaths
#'
#' Returns a vector containing all the file paths for positive samples.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return pos_file_paths. A vector containing file paths for positive samples.
getAllPosGSKFilePaths <- function() {
    pos_files <- meta_all[, "file_name_pos"]
    pos_files <- as.character(pos_files)
    pos_file_paths <- character()
    for (i in 1:length(pos_files)) {
        if (grepl("block1", pos_files[i])) {
            neg_file_path <- paste(neg_dir, "/block1pos/", pos_files[i], ".cdf", sep = "")
            pos_file_paths[i] <- neg_file_path
        }
        else if (grepl("block2", pos_files[i])) {
            neg_file_path <- paste(neg_dir, "/block2pos/", pos_files[i], ".cdf", sep = "")
            pos_file_paths[i] <- neg_file_path
        }
        else if (grepl("block3", pos_files[i])) {
            neg_file_path <- paste(neg_dir, "/block3pos/", pos_files[i], ".cdf", sep = "")
            pos_file_paths[i] <- neg_file_path
        }
        else if (grepl("block4", pos_files[i])) {
            neg_file_path <- paste(neg_dir, "/block4pos/", pos_files[i], ".cdf", sep = "")
            pos_file_paths[i] <- neg_file_path
        }
        else {
            pos_file_paths[i] <- ""
        }
    }
    return(pos_file_paths)
}

#' getNegGSKFilePath
#'
#' Returns the file path for a given name of a sample measured in negative mode.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return gsk_neg_file_path. A vector containing block metadata.
getNegGSKFilePath <- function(gsk_filename) {
    if (grepl("block1", gsk_filename)) {
        gsk_neg_file_path <- paste(neg_dir, "/block1neg/", gsk_filename, ".cdf", sep = "")
    }
    else if (grepl("block2", gsk_filename)) {
        gsk_neg_file_path <- paste(neg_dir, "/block2neg/", gsk_filename, ".cdf", sep = "")
    }
    else if (grepl("block3", gsk_filename)) {
        gsk_neg_file_path <- paste(neg_dir, "/block3neg/", gsk_filename, ".cdf", sep = "")
    }
    else if (grepl("block4", gsk_filename)) {
        gsk_neg_file_path <- paste(neg_dir, "/block4neg/", gsk_filename, ".cdf", sep = "")
    }
    else {
        gsk_neg_file_path <- ""
    }
    return(gsk_neg_file_path)
}

#' getPosGSKFilePath
#'
#' Returns the file path for a given name of a sample measured in positive mode.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return gsk_pos_file_path. A vector containing block metadata.
getPosGSKFilePath <- function(gsk_filename) {
    if (grepl("block1", gsk_filename)) {
        gsk_pos_file_path <- paste(pos_dir, "/block1pos/", gsk_filename, ".cdf", sep = "")
    }
    else if (grepl("block2", gsk_filename)) {
        gsk_pos_file_path <- paste(pos_dir, "/block2pos/", gsk_filename, ".cdf", sep = "")
    }
    else if (grepl("block3", gsk_filename)) {
        gsk_pos_file_path <- paste(pos_dir, "/block3pos/", gsk_filename, ".cdf", sep = "")
    }
    else if (grepl("block4", gsk_filename)) {
        gsk_pos_file_path <- paste(pos_dir, "/block4pos/", gsk_filename, ".cdf", sep = "")
    }
    else {
        gsk_pos_file_path <- ""
    }
    return(gsk_pos_file_path)
}

#' getBlockMetadata
#'
#' Returns the block information assoicated with a vector
#' of sample names.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return block_metadata. A vector containing block metadata.
getBlockMetadata <- function(sample_names) {
    block_metadata <- character()
    for (i in 1:length(sample_names)) {
        if (grepl("block1", sample_names[i])) {
            block_metadata[i] <- "block1"
        }
        else if (grepl("block2", sample_names[i])) {
            block_metadata[i] <- "block2"
        }
        else if (grepl("block3", sample_names[i])) {
            block_metadata[i] <- "block3"
        }
        else if (grepl("block4", sample_names[i])) {
            block_metadata[i] <- "block4"
        }
        else {
            block_metadata[i] <- ""
        }
    }
    return(block_metadata)
}

#' getQCSampleMetadata
#'
#' Returns the sample type information assoicated with a vector
#' of sample names.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return sample_meta. A vector containing sample type metadata.
getQCSampleMetadata <- function(samples, esi_mode="negative") {
    col_name="file_name_neg"
    if (esi_mode == "positive") {
        col_name == "file_name_pos"
    }
    sample_meta <- character(0)
    for (i in 1:length(samples)) {
        if (meta_all[which(meta_all[, col_name]==samples[i]), "type"] == 'QC') {
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
getRegimenMetadata <- function(samples, esi_mode="negative") {
    col_name="file_name_neg"
    if (esi_mode == "positive") {
        col_name == "file_name_pos"
    }

    sample_meta <- character(0)
    for (i in 1:length(samples)) {
        regimen <- meta_all[which(meta_all[, col_name]==samples[i]), "Regimen"]
        if (is.na(regimen)) {
            sample_meta[i] <- "QC"
        }
        else {
            sample_meta[i] <- regimen
        }
    }
    return(sample_meta)
}

#' getNegativeQCSampleNames
#'
#' Returns all QC sample names which were analysed in negative ionization
#' mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return neg_qc_names. A vector containing negative QC sample names.
getNegativeQCSampleNames <- function() {
    neg_qc_names <- meta_all[meta_all[, "type"] == 'QC', "file_name_neg"]
    neg_qc_names <- as.character(neg_qc_names)
    return(neg_qc_names)
}

#' getPositiveQCSampleNames
#'
#' Returns all QC sample names which were analysed in positive ionization
#' mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return pos_qc_names. A vector containing positive QC sample names.
getPositiveQCSampleNames <- function() {
    pos_qc_names <- meta_all[meta_all[, "type"] == 'QC', "file_name_pos"]
    pos_qc_names <- as.character(pos_qc_names)
    return(pos_qc_names)
}

#' getNumberOfNegativeQCs
#'
#' Returns the total number of QCs measured in negative ionization mode
#  in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer representing the number negative QCs.
getNumberOfNegativeQCs <- function() {
    neg_qc_names <- getNegativeQCSampleNames()
    return(length(neg_qc_names))
}

#' getNumberOfPositiveQCs
#'
#' Returns the total number of QCs measured in positive ionization mode
#  in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer represe1nting the number positive QCs.
getNumberOfPositiveQCs <- function() {
    pos_qc_names <- getPositiveQCSampleNames()
    return(length(pos_qc_names))
}

#' getNegativePlasmaSampleNames
#'
#' Returns a vector of plasma sample names which have been measured in
#' negative ionization mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return neg_plasma_names. A vector containing a list of negative plasma sample names.
getNegativePlasmaSampleNames <- function() {
    neg_plasma_names <- meta_all[meta_all[, "type"] == 'Sample', "file_name_neg"]
    neg_plasma_names <- as.character(neg_plasma_names)
    return(neg_plasma_names)
}

#' getPositivePlasmaSampleNames
#'
#' Returns a vector of plasma sample names which have been measured in
#' positive ionization mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return pos_plasma_names. A vector containing a list of positive plasma sample names.
getPositivePlasmaSampleNames <- function() {
    pos_plasma_names <- meta_all[meta_all[, "type"] == 'Sample', "file_name_pos"]
    pos_plasma_names <- as.character(pos_plasma_names)
    return(pos_plasma_names)
}

#' getNumberOfNegativePlasmaSamples
#'
#' Returns an integer showing the number of plasma sample names which
#' have been measured in negative ionization mode in the GSK
#' longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer showing number of negative plasma sample names.
getNumberOfNegativePlasmaSamples <- function() {
    neg_plasma_names <- getNegativePlasmaSampleNames()
    return(length(neg_plasma_names))
}

#' getNumberOfPositivePlasmaSamples
#'
#' Returns an integer showing the number of plasma sample names which
#' have been measured in positive ionization mode in the GSK
#' longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer showing number of positive plasma sample names.
getNumberOfPositivePlasmaSamples <- function() {
    pos_plasma_names <- getPositivePlasmaSampleNames()
    return(length(pos_plasma_names))
}

#' getNegativeTP0PlasmaSampleNames
#'
#' Returns an integer showing the number of plasma sample names which
#' have been measured at Time Point 0 in negative ionization mode in
#' the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return tp0_sample_names. An integer showing number of positive plasma sample names.
getNegativeTP0PlasmaSampleNames <- function() {
    tp0_sample_names <- meta_all[meta_all[, "Timepoint..PIMS."] == 0, "file_name_neg"]
    tp0_sample_names <- as.character(tp0_sample_names)
    tp0_sample_names <- na.omit(tp0_sample_names)
    return(tp0_sample_names)
}

#' getPositiveTP0PlasmaSampleNames
#'
#' Returns an integer showing the number of plasma sample names which
#' have been measured at Time Point 0 in negative ionization mode in
#' the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return tp0_sample_names. An integer showing number of positive plasma sample names.
getPositiveTP0PlasmaSampleNames <- function() {
    tp0_sample_names <- meta_all[meta_all[, "Timepoint..PIMS."] == 0, "file_name_pos"]
    tp0_sample_names <- as.character(tp0_sample_names)
    tp0_sample_names <- na.omit(tp0_sample_names)
    return(tp0_sample_names)
}