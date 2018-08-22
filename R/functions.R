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

#' getAllGSKFilePaths
#'
#' Returns a vector containing all the file paths for negative or positive samples.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return file_paths. A vector containing file paths for negative or positive samples.
getAllGSKFilePaths <- function(mode="negative") {
    file_paths <- character()

    if (mode == "negative") {
        files <- meta_all[, "file_name_neg"]
        files <- as.character(files)
        for (i in 1:length(files)) {
            if (grepl("block1", files[i])) {
                file_paths[i] <- paste(neg_dir, "/block1neg/", files[i], ".cdf", sep = "")
            }
            else if (grepl("block2", files[i])) {
                file_paths[i] <- paste(neg_dir, "/block2neg/", files[i], ".cdf", sep = "")
            }
            else if (grepl("block3", files[i])) {
                file_paths[i] <- paste(neg_dir, "/block3neg/", files[i], ".cdf", sep = "")
            }
            else if (grepl("block4", files[i])) {
                file_paths[i] <- paste(neg_dir, "/block4neg/", files[i], ".cdf", sep = "")
            }
        }
    }
    else {
        files <- meta_all[, "file_name_pos"]
        files <- as.character(files)
        for (i in 1:length(files)) {
            if (grepl("block1", files[i])) {
                file_paths[i] <- paste(pos_dir, "/block1pos/", files[i], ".cdf", sep = "")
            }
            else if (grepl("block2", files[i])) {
                file_paths[i] <- paste(pos_dir, "/block2pos/", files[i], ".cdf", sep = "")
            }
            else if (grepl("block3", files[i])) {
                file_paths[i] <- paste(pos_dir, "/block3pos/", files[i], ".cdf", sep = "")
            }
            else if (grepl("block4", files[i])) {
                file_paths[i] <- paste(pos_dir, "/block4pos/", files[i], ".cdf", sep = "")
            }
        }
    }
    return(file_paths)
}

#' getGSKFilePath
#'
#' Returns the file path for a given name of a sample measured in negative or positive mode.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return gsk_file_path. A vector containing block metadata.
getGSKFilePath <- function(gsk_filename, mode="negative") {
    gsk_file_path <- ""
    if (mode == "negative") {
        if (grepl("block1", gsk_filename)) {
            gsk_file_path <- paste(neg_dir, "/block1neg/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block2", gsk_filename)) {
            gsk_file_path <- paste(neg_dir, "/block2neg/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block3", gsk_filename)) {
            gsk_file_path <- paste(neg_dir, "/block3neg/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block4", gsk_filename)) {
            gsk_file_path <- paste(neg_dir, "/block4neg/", gsk_filename, ".cdf", sep = "")
        }
    }
    else {
        if (grepl("block1", gsk_filename)) {
            gsk_file_path <- paste(pos_dir, "/block1pos/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block2", gsk_filename)) {
            gsk_file_path <- paste(pos_dir, "/block2pos/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block3", gsk_filename)) {
            gsk_file_path <- paste(pos_dir, "/block3pos/", gsk_filename, ".cdf", sep = "")
        }
        else if (grepl("block4", gsk_filename)) {
            gsk_file_path <- paste(pos_dir, "/block4pos/", gsk_filename, ".cdf", sep = "")
        }
    }
    return(gsk_file_path)
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

#' getQCSampleNames
#'
#' Returns all QC sample names which were analysed in negative or positive ionization
#' mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return qc_names. A vector containing negative or positive QC sample names.
getQCSampleNames <- function(mode="negative") {
    if (mode == "negative") {
        col_name <- "file_name_neg"
    }
    else {
        col_name <- "file_name_pos"
    }

    qc_names <- meta_all[meta_all[, "type"] == 'QC', col_name]
    qc_names <- as.character(qc_names)
    return(qc_names)
}

#' getNumberOfQCs
#'
#' Returns the total number of QCs measured in negative or positive ionization mode
#  in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer representing the number negative or positive QCs.
getNumberOfQCs <- function(mode="negative") {
    qc_names <- getQCSampleNames(mode)
    return(length(qc_names))
}

#' getPlasmaSampleNames
#'
#' Returns a vector of plasma sample names which have been measured in
#' negative or positive ionization mode in the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return plasma_names. A vector containing a list of negative plasma sample names.
getPlasmaSampleNames <- function(mode="negative") {
    if (mode == "negative") {
        col_name <- "file_name_neg"
    }
    else {
        col_name <- "file_name_pos"
    }
    plasma_names <- meta_all[meta_all[, "type"] == 'Sample', col_name]
    plasma_names <- as.character(plasma_names)
    return(plasma_names)
}

#' getNumberOfPlasmaSamples
#'
#' Returns an integer showing the number of plasma sample names which
#' have been measured in negative or positive ionization mode in the GSK
#' longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return An integer showing number of negative or positive plasma sample names.
getNumberOfPlasmaSamples <- function(mode="negative") {
    plasma_names <- getPlasmaSampleNames(mode)
    return(length(plasma_names))
}

#' getTP0PlasmaSampleNames
#'
#' Returns an integer showing the number of plasma sample names which
#' have been measured at Time Point 0 in negative or positive ionization mode in
#' the GSK longitudinal data set.
#' @author Peter Li <peter@gigasciencejournal.com>
#' @return tp0_sample_names. An integer showing number of positive or negative plasma sample names.
getNegativeTP0PlasmaSampleNames <- function(mode="negative") {
    if (mode == "negative") {
        col_name <- "file_name_neg"
    }
    else {
        col_name <- "file_name_pos"
    }
    tp0_sample_names <- meta_all[meta_all[, "Timepoint..PIMS."] == 0, col_name]
    tp0_sample_names <- as.character(tp0_sample_names)
    tp0_sample_names <- na.omit(tp0_sample_names)
    return(tp0_sample_names)
}
