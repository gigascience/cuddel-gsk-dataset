# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 14/8/2018

# Location of GSK data set
datadir = "/home/peter/"
# File path for negative files
neg_dir = paste(datadir, "gsk/raw/esi_neg/netcdf", sep ="")

######################################
# Read in metadata for data analysis #
######################################

# Read in and sort metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta.csv", sep =""))
meta <- meta[naturalorder(meta$file_name_neg),]
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all.csv", sep =""))
meta_all <- meta_all[naturalorder(meta_all$file_name_neg),]

getAllNegGSKFilePaths <- function() {
    neg_files <- meta_all[,"file_name_neg"]
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


getQCSampleMetadata <- function(samples) {
    sample_meta <- character(0)
    for (i in 1:length(samples)) {
        if (meta_all[which(meta_all[,"file_name_neg"]==samples[i]), "type"] == 'QC') {
            sample_meta[i] <- "QC"
        }
        else {
            sample_meta[i] <- "Sample"
        }
    }
    return(sample_meta)
}