# Title     : metab_anno.R
# Objective : Annotate peaks with metabolite identities using MAIT
# Created by: peterli
# Created on: 30/10/2018

library(xcms)
library(CAMERA)
library(MAIT)

source("functions.R")

# cd to /home/peter/gsk/raw/esi_pos/netcdf

# Load positive data
setwd("/home/peter/gsk/raw/esi_pos/netcdf")

# These commands were executed using XCMS and CAMERA
######################################
# Apply XCMS onto QC and sample data #
######################################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and takes time to complete!!
pos_xset <- xcmsSet(pos_file_paths, step=0.02, snthresh=3, mzdiff=0.05)
# Match peaks representing same analyte across samples
grp_pos_xset <- group(pos_xset, bw=10, mzwid=0.05)
# Create CAMERA object
pos_xsa <- xsAnnotate(grp_pos_xset)
# Get peaklist - contains 6200 peaks identified by XCMS
pos_peaklist <- getPeaklist(pos_xsa)

# Get peaks, masses and rt data
sample_names <- grep("GSK_pos", names(pos_peaklist), value=TRUE)
peaks <- pos_peaklist[, sample_names]
peaks <- as.matrix(peaks)
masses <- getPeaklist(pos_xsa)$mz
rt <- getPeaklist(pos_xsa)$rt/60

# Perform an annotation stage and metabolite identification
importMAIT <- MAITbuilder(data=peaks, masses=masses,
rt=rt, significantFeatures=TRUE,
spectraEstimation=TRUE, rtRange=0.2,
corThresh=0.7)

# Improve peak annotation using biotransformations function with positive mode data
importMAIT <- Biotransformations(MAIT.object=importMAIT,
adductAnnotation=TRUE, peakPrecision=0.005, adductTable=NULL)

# Identify metabolites
importMAIT <- identifyMetabolites(MAIT.object=importMAIT,
peakTolerance=0.005, polarity="positive")

# Read in metabolite annotations
mets <- read.csv("Tables/metaboliteTable.csv")
head(mets)