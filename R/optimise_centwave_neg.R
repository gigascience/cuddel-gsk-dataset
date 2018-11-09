# Title     : optimise_centwave_neg.R
# Objective : TODO
# Created by: peterli
# Created on: 7/11/2018

# optimise_centwave_neg.R
#
# Based on https://github.com/vanmooylipidomics/LipidomicsToolbox/blob/master/prepOrbidata.R
#
# Created 11/18/2015 by Jamie Collins, james.r.collins@aya.yale.edu
# Name changed from LOBSTAHS_main.R to prepOrbidata.R by J.R.C., 12/4/2015
#
# Purpose: Use xcms, CAMERA, and (optionally) IPO to prepare a multiple sample,
# HPLC-ESI-MS dataset from the Exactive Plus Orbitrap for follow-on feature ID
# and annotation with the LOBSTAHS (Lipid and Oxylipin Biomarker Screening
# through Adduct Hierarchy Sequences) package.
#
# LOBSTAHS is under development by the Van Mooy Lab at Woods Hole Oceanographic
# Institution. Currently, this script is written to analyze lipid data from the
# experiment described in Graff van Creveld et al., 2015, "Early perturbation in
# mitochondria redox homeostasis in response to environmental stress predicts
# cell fate in diatoms," ISME Journal 9:385-395. This dataset is used to
# demonstrate the use of LOBSTAHS in Collins, J.R., B.R. Edwards, H.F.
# Fredricks, and B.A.S. Van Mooy, "Semi-untargeted discovery and identification
# of oxidative stress biomarkers using a lipidomics pipeline for complex
# datasets."
#
# This script:
#
# 1. Uses xcms to perform (1) peak picking and integration, (2) chromatographic
# alignment, and (3) nonlinear grouping across samples. Requires package "IPO"
# for parameter optimization.
#
# 2. Uses CAMERA to perform (1) identification of secondary isotope peaks,
# (2) creation of CAMERA pseudospectra using correlation xcms peak groups
# between and within samples, and (3) creation of a CAMERA xsAnnotate object
# suitable for input to the LOBSTAHS function "doLOBscreen"
#
# As described in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van
# Mooy, "Untargeted discovery and identification of oxidative stress biomarkers
# using a lipidomics pipeline for complex datasets"
#
# Obtain current versions of packages, scripts, and necessary dependencies at
# https://github.com/vanmooylipidomics/
#
# Please direct questions/suggestions/comments to Jamie Collins,
# james.r.collins@aya.yale.edu, or Helen Fredricks, hfredricks@whoi.edu
#
# Revision history maintained on GitHub
#
################ Caveats and prerequisites #############
#
# Presumes user has installed the R packages "xcms", "CAMERA", "tools", "IPO",
# with all required dependencies
#
# If multicore tasking is desired, "snowfall" also required; Rmpi doesn't seem
# to be necessary
#
# This script the following inputs:
#
# 1. A series of .mzXML files from the same dataset, containing centroided ms1
# data of a single ion mode. File conversion from the Thermo .raw format,
# centroiding of data, and extraction of + and - mode scans into separate files
# can be accomplished in batch using the script
# "Exactive_full_scan_process_ms1+.r", available from
# https://github.com/vanmooylipidomics/LipidomicsToolbox. The .mzXML files
# should be placed together in a single directory, which can be specified by the
# user below.
#
# 2. If the package IPO was previously used to optimize xcms peak-picking or
# group/retcor parameters AND automatic import of the optimized settings from an
# existing .csv file is desired, specification of the path to file
# "IPO_xcmsparamfits_ ... .csv," where ... is an ISO 8601 timestamp. A suitable
# .csv file will be generated if the user elects IPO at two user-input points in
# this script, or such a file can be generated from IPO using the helper script
# optim_centWaveParams_standalone.R, latest version at
# https://github.com/vanmooylipidomics/LipidomicsToolbox/blob/master/optim_centWaveParams_standalone.R

# Initial setup and variable definition

# Load required packages
library(tools)
library(xcms)
library(CAMERA)
library(rsm)
library(IPO)
library(snow)  # For multicore tasking

# Load helper functions for selecting data files
source('functions.R')

# ******************************************************************
################ Basic user begin editing here #############
# ******************************************************************

# Use functions.R to specify negative QC files for IPO to use
neg_qc_sample_names <- getQCSampleNames()
neg_qc_files <- lapply(neg_qc_sample_names, getGSKFilePath)
neg_qc_files <- unlist(neg_qc_files)

# User: define locations of data files
wd <- "/home/peter/gsk/raw/esi_neg/netcdf/"
setwd(wd)

# Specify directories subordinate to the working directory in which the .mzXML
# files for xcms can be found; per xcms documentation, use subdirectories within
# these to divide files according to treatment/primary environmental variable
# (e.g., station number along a cruise transect) and file names to indicate
# timepoint/secondary environmental variable (e.g., depth)
block_dirs <- c("block1neg/", "block2neg/", "block3neg", "block4neg")

# Specify which of the directories above you wish to analyze this time through
chosenFileSubset = "Pt_H2O2_mzXML_ms1_pos/"

# Specify files you don't want to push through xcms (e.g., blanks)
excluded_files = c("0475", "0474")

# If planning to use IPO, specify the ID numbers (i.e., Orbi_xxxx.mzXML) of the
# files you'd like to use for optimization; otherwise, IPO will try to use the
# entire dataset and you'll probably wind up with errors
IPO.filesubset = c("0468", "0476", "0477")

# If you aren't planning on running IPO to optimize centWave and/or group/retcor
# parameters this session, but you have some parameter values from an earlier
# IPO run saved in a .csv file, you can specify the file paths below. you will
# still be given the option later to choose which parameters to use.
saved_IPO_params_centW = "IPO_xcmsparamfits_2018-11-09T18-36-10-0300.csv"
saved_IPO_params_groupretcor = "IPO_retcorGroupparamfits_2018-11-09T18-36-10-0300.csv"

# Specify which xcms centWave parameter values to use, and whether to run IPO

# Run IPO for optimization now, then use those settings
centWparam.source <- 1
# Use default settings specified in script below
## centWparam.source = 2
# Read in previously optimized parameter values from the .csv file specified
# immediately above as saved_IPO_params_centW
## centWparam.source = 3

# user: specify which xcms group and retcor parameter values to use, and whether to run IPO

# groupretcor.prams.source = 1 # to run IPO for optimization now, then use those settings
groupretcor.prams.source = 2 # to use default settings specified in script below
# groupretcor.prams.source = 3 # to read in previously optimized parameter values from the .csv file specified immediately above as saved_IPO_params_groupretcor

# Specify which retcor method to use: 'loess' or 'obiwarp'
retcor.meth = "loess"
# retcor.meth = "obiwarp"

# ******************************************************************
################ Basic user stop editing here #############
# ******************************************************************

################# Define functions; run me first #############

# readinteger: for a given prompt, allows capture of user input as an integer;
# rejects non-integer input
readinteger = function(prompttext) {
    n = readline(prompt=prompttext)
    if (!grepl("^[0-9]+$", n)) {
        return(readinteger(prompttext))
    }
    as.integer(n)
}

# readyesno: for a given prompt, allows capture of user input as y or n; rejects
# other input
readyesno = function(prompttext) {
    n = readline(prompt=prompttext)
    if (!grepl("y|n", n)) {
        return(readyesno(prompttext))
    }
    as.character(n)
}

# verifyFileIonMode: return the ion mode of data in a particular mzXML file, by
# examining "polarity" attribute of each scan in the file
verifyFileIonMode = function(mzXMLfile) {
    rawfile = xcmsRaw(mzXMLfile) # create an xcmsraw object out of the first file
    # determine ion mode by examining identifier attached to scan events
    if (table(rawfile@polarity)["negative"]==0 & (table(rawfile@polarity)["positive"]==length(rawfile@scanindex))) { # can deduce that the file contains positive mode data
        filepolarity = 1 # positive
    } else if (table(rawfile@polarity)["positive"]==0 & (table(rawfile@polarity)["negative"]==length(rawfile@scanindex))) { # probably negative mode data
        filepolarity = -1 # negative
    } else if (table(rawfile@polarity)["positive"]>=1 & table(rawfile@polarity)["negative"]>=1) { # scans of both mode present in the file; the original .raw files weren't split by mode during initial .mzXML conversion, or something else is wrong
        stop("At least one file in the current dataset contains scans of more than one ion mode. Please ensure data for different ion modes have been extracted into separate files. Stopping...") # stop script if this is the case
    } else if (table(rawfile@polarity)["positive"]==0 & table(rawfile@polarity)["negative"]==0) {
        stop("Can't determine ion mode of data in the first file. Check manner in which files were converted. Stopping...") # stop script if this is the case
    }
    filepolarity
}

# getSubsetIonMode: return the ion mode of a subset of files, using sapply of
# verifyFileIonMode
getSubsetIonMode = function(mzXMLfilelist) {
    ionmodecount = sum(sapply(mzXMLfilelist, verifyFileIonMode)) # get sum of ion mode indicators for the files in the subset
    if (ionmodecount==length(mzXMLfilelist)) { # can conclude that all files contain positive mode data
        subset.polarity = "positive"
    } else if (ionmodecount==-length(mzXMLfilelist)) { # can conclude that all files contain negative mode data
        subset.polarity = "negative"
    }
    subset.polarity
}

# selectXMLSubDir: allows user to choose which subset of files to process
selectXMLSubDir = function(mzXMLdirList) {
    print(paste0("mzXML files exist in the following directories:"))
    for (i in 1:length(mzXMLdirList)) {
        # get number of mzXML files in this directory
        numGoodFiles = length(list.files(mzXMLdirList[i],
        recursive = TRUE, full.names = TRUE,
        pattern = "*(.mzXML|.mzxml)"))

        if (numGoodFiles>0) { # there are .mzXML data files in this directory
            print(paste0(i, ". ", numGoodFiles," .mzXML files in directory '",mzXMLdirList[i],"'"))
        }
    }

    processDecision = readinteger("Specify which subset you'd like to process, using integer input: ")
    mzXMLdirList[processDecision]
}

# getFNmatches: returns index(es) of file names in a given file list containing
# the ID numbers in a match list
getFNmatches = function(filelist,IDnumlist) {
    unique(grep(paste(IDnumlist,collapse="|"),filelist, value=FALSE))
}

# genTimeStamp: generates a timestamp string based on the current system time
genTimeStamp = function () {
    output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
    output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
    output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)
}

# check to make sure user has specified at least something in block_dirs
if (!exists("block_dirs")) {
    stop("User has not specified any directories containing mzXML files. Specify a value for block_dirs.")
}

################# Load in mzXML files, get xcms settings from IPO or user input #############
# load selected subset for processing
mzXMLfiles.raw = list.files(chosenFileSubset, recursive = TRUE, full.names = TRUE)
# verify the ion mode of the data in these files
subset.polarity = getSubsetIonMode(mzXMLfiles.raw)
# provide some feedback to user
print(paste0("Loaded ",
    length(mzXMLfiles.raw),
    " mzXML files. These files contain ",
    subset.polarity," ion mode data. Raw dataset consists of:"))

print(mzXMLfiles.raw)

# check whether user has elected to exclude any files, and exclude them if they happen to be in this subset
if (exists("excluded_files") & length("excluded_files")>0) {
    # index files to be excluded
    excludedfiles <- getFNmatches(IDnumlist = excluded_files, filelist = mzXMLfiles.raw)
    print(paste0("The following files will be excluded from processing based on user's input:"))
    print(mzXMLfiles.raw[excludedfiles])
    # exclude the files from mzXMLfiles
    mzXMLfiles = mzXMLfiles.raw[-excludedfiles]
} else {
    mzXMLfiles = mzXMLfiles.raw
}


# Peak-picking & creation of xcmsSet using xcms (and IPO, if desired)
if (centWparam.source==1) {  # Use IPO to optimize settings for centWave method
    # Define ranges of parameters for testing. If a single value is specified
    # or centWave default is used, parameter will not be optimized
    peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')  # get defaults
    peakpickingParameters$min_peakwidth <- c(7.5) # centerpoint is 15
    peakpickingParameters$max_peakwidth <- c(26) # centerpoint is 60
    # Set ppm low to avoid peak data insertion errors from centWave; IPO wants
    # to use something like 5.5 ppm if you allow it to "optimize," but this is
    # too high. Jan Stanstrup says ppm = 20 for orbitrap data, see:
    # http://www.metabolomics-forum.com/index.php?topic=1246.0
    peakpickingParameters$ppm <- c(20)  # Was originally 2.5 in this script
    # A very long optimization routine settled on a value of 2.4 for prefilter
    peakpickingParameters$prefilter <- 3
    peakpickingParameters$value_of_prefilter <- c(100, 5000)
    peakpickingParameters$snthresh <- c(10)
    peakpickingParameters$noise <- c(500)

    print(paste0("Using R package IPO to optimize centWave peak-picking settings with starting parameters user has specified in script. Using following subset of files for optimization:"))

    if (exists("IPO.filesubset") & length("IPO.filesubset")>0) {
        # Get indexes to subset
        IPOsubset <- mzXMLfiles[getFNmatches(IDnumlist = IPO.filesubset, filelist = mzXMLfiles)]
        print(IPOsubset)
    } else {
        print(paste0("User did not specify a subset of files for IPO optimization."))
        print(paste0("Defaulting to all files in the current dataset. This will take long time and may yield errors."))
        IPOsubset <- mzXMLfiles
    }

    resultPeakpicking <- optimizeXcmsSet(files=neg_qc_files[1:5],
        params=peakpickingParameters,
        BPPARAM = MulticoreParam(workers = detectCores()),
        subdir=NULL)

    optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset

    #### Export IPO starting value(s) and optimal settings for each parameter to .csv ####

    # Write 3-column table to .csv using write.table()
    peakpicking.exportmat <- cbind(sort(rownames(as.matrix(peakpickingParameters))),
        as.character(resultPeakpicking$best_settings$parameters[sort(names(resultPeakpicking$best_settings$parameters))]),
        as.character(peakpickingParameters[sort(rownames(as.matrix(peakpickingParameters)))]))

    # Append an additional row with the file names of the files used for
    # optimization
    peakpicking.exportmat <- rbind(peakpicking.exportmat,
        c("Files_used_for_optimization",
        paste(IPOsubset, collapse = ", "),
        ""))

    timestamp.now <- print(genTimeStamp())

    write.table(peakpicking.exportmat,
        file = paste("IPO_centWaveparamfits_", timestamp.now, ".csv", sep=""),
        col.names = c("centWave_parameter", "IPO_optim_value", "Starting_value(s)"),
        row.names = FALSE,
        sep=",")

    print(paste0("IPO optimization complete. Optimized and starting values for"))
    print(paste0("centWave parameters written to file: IPO_centWaveparamfits_", timestamp.now,".csv"))

    # Use the just-obtained optimized parameter values for xcmsSet creation
    print(paste0("Using IPO-optimized settings for findPeaks.centWave..."))
    centW.min_peakwidth <- resultPeakpicking$best_settings$parameters$min_peakwidth
    centW.max_peakwidth <- resultPeakpicking$best_settings$parameters$max_peakwidth
    centW.ppm           <- resultPeakpicking$best_settings$parameters$ppm
    centW.mzdiff        <- resultPeakpicking$best_settings$parameters$mzdiff
    centW.snthresh      <- resultPeakpicking$best_settings$parameters$snthresh
    centW.prefilter     <- c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter)
    centW.noise         <- resultPeakpicking$best_settings$parameters$noise

    # Not using IPO settings from resultPeakpicking$best_settings$parameters for
    # mzCenterFun, integrate, fitgauss, verbose.columns, nSlaves since those
    # weren't targets of optimization

} else if (centWparam.source==2) { # user wants to use settings specified below

    print(paste0("Using values of centWave parameters specified in the script by user..."))

    # "non-optimized" settings listed here are based on recommended
    # "HPLC/Orbitrap settings" from Table 1 of Patti et al., 2012,
    # "Meta-analysis of untargeted metabolomic data from multiple profiling
    # experiment," Nature Protocols 7: 508-516

    centW.min_peakwidth <- 10
    # lowered from Patti et al. recommended HPLC setting of 60 based on visual
    # inspection of a single sample with plotPeaks
    centW.max_peakwidth <- 45
    centW.ppm <- 2.5
    centW.mzdiff <- 0.005
    centW.snthresh <- 10
    centW.prefilter <- c(3,7500) # 3.5k recommended by Patti et al. appears to be too low
    centW.noise <- 500

} else if (centWparam.source==3) { # user wants to read in parameter values from file

    print(paste0("Loading values of centWave parameters from previous IPO optimization run in .csv file:"))
    print(paste0(saved_IPO_params_centW))

    centWprams.from.file = read.csv(saved_IPO_params_centW,colClasses = "character")

    centW.min_peakwidth <- as.numeric(centWprams.from.file[centWprams.from.file[,1]=="min_peakwidth",2])
    centW.max_peakwidth <- as.numeric(centWprams.from.file[centWprams.from.file[,1]=="max_peakwidth",2])
    centW.ppm           <- as.numeric(centWprams.from.file[centWprams.from.file[,1]=="ppm",2])
    centW.mzdiff        <- as.numeric(centWprams.from.file[centWprams.from.file[,1]=="mzdiff",2])
    centW.snthresh      <- as.numeric(centWprams.from.file[centWprams.from.file[,1]=="snthresh",2])
    centW.prefilter     <- c(as.numeric(centWprams.from.file[centWprams.from.file[,1]=="prefilter",2]),as.numeric(centWprams.from.file[centWprams.from.file[,1]=="value_of_prefilter",2]))
    centW.noise         <- as.numeric(centWprams.from.file[centWprams.from.file[,1]=="noise",2])

}

# specify some additional settings we wish to keep constant, regardless of where
# the parameters above were obtained
centW.fitgauss = TRUE
centW.sleep = 1
centW.mzCenterFun = c("wMean")
centW.verbose.columns = TRUE
centW.integrate = 1
# setting this very low, per Jan Stanstrup; low setting uses more memory but
# helps avoid the situation where mass accuracy eclipses the actual width of the
# m/z windows used to define each peak (a real possibility with Orbitrap data;
# see http://metabolomics-forum.com/viewtopic.php?f=8&t=598#p1853)
centW.profparam = list(step=0.001)
# if you have r package "snow" installed, can set to number of cores you wish to
# make use of
centW.nSlaves = 4

# ################# Peak visualization using individual sample files #############

# # optional section for method development

# # create xcmsRaw object from just a single sample (for method development)

# xfile_raw = xcmsRaw(mzXMLfiles[1], profparam = centW.profparam, includeMSn=TRUE)
# profStep(xfile_raw) = 0.005

# rawpeaks = findPeaks.centWave(xfile_raw,
# ppm = centW.ppm,
# peakwidth = c(centW.min_peakwidth,centW.max_peakwidth),
# fitgauss = centW.fitgauss,
# noise = centW.noise,
# mzdiff = centW.mzdiff,
# verbose.columns = centW.verbose.columns,
# snthresh = centW.snthresh,
# integrate = centW.integrate,
# prefilter = centW.prefilter,
# mzCenterFun = centW.mzCenterFun
# #                 ,sleep = centW.sleep
# #                 nSlaves = centW.nSlaves
# )

# # # despite the good press, massifquant was picking some very bad looking features, using centWave for time being instead

# # # rawpeaks = findPeaks.massifquant(xfile_raw,
# # criticalValue = 1,
# # consecMissedLimit = 2, # supposedly optimal for Orbitrap data
# # prefilter = c(3,10000), # documentation says the first argument only necessary if withWave = 1, but then the example shows it there with withWave = 0; using 10k rather than the 3.5k recommended by Patti et al.
# # ppm = 2.5, # using recommended setting of Patti et al., 2012 (see below)
# # unions = 1,
# # profparam = centW.profparam,
# # withWave = 1, # two arguments immediately below are if withWave = 1 only
# # sleep = 1,
# # peakwidth = c(10,45), # min. feature length in time scans, max chromatographic peak width
# # snthresh = 10,
# # integrate = 1,
# # checkBack = 1,
# # fitgauss = FALSE
# # #                nSlaves = 4 # if you have r package "snow" installed, can set to number of cores you wish to make use of
# # )

# # plot some selected peaks

# plotPeaks(xfile_raw,rawpeaks[10150:10174,],figs = c(5,5),width = 100)
# plotPeaks(xfile_raw,rawpeaks[150:174,],figs = c(5,5),width = 100)
# plotPeaks(xfile_raw,rawpeaks[1:24,],figs = c(5,5),width = 100)

# # N.B., just because you can't see the full extent of the peaks in some of the subplots doesn't mean they're bad; appears to be something wonky with the ylim setting in plotPeaks; see www.metabolomics-forum.com/viewtopic.php?f=8&t=875

# # for example, can look at an individual peak this way:

# plotEIC(xfile_raw, mzrange = rawpeaks[10150,c("mzmin","mzmax")], rtrange = rawpeaks[10150,c("rtmin","rtmax")]   )

# # a diagnostic plot, showing necessity of setting profparam low enough

# mz_width = rawpeaks@.Data[,"mzmax"] - rawpeaks@.Data[,"mzmin"]
# plot(density(mz_width,adjust=0.2))

# # some other plots

# plotChrom(xfile_raw)

################# Create xcmsSet using selected settings #############

print(paste0("Creating xcmsSet object from ",length(mzXMLfiles)," mzXML files remaining in dataset using specified settings..."))

# create xcms xset object; runs WAY faster with multicore tasking enabled;

xset_centWave = xcmsSet(mzXMLfiles,
method = "centWave",
profparam = centW.profparam,
ppm = centW.ppm,
peakwidth = c(centW.min_peakwidth,centW.max_peakwidth),
fitgauss = centW.fitgauss,
noise = centW.noise,
mzdiff = centW.mzdiff,
verbose.columns = centW.verbose.columns,
snthresh = centW.snthresh,
integrate = centW.integrate,
prefilter = centW.prefilter,
mzCenterFun = centW.mzCenterFun,
#                 sleep = centW.sleep
nSlaves = centW.nSlaves)

print(paste0("xcmsSet object xset_centWave created:"))

print(xset_centWave)

# Some notes:
#
#  1. If using massifquant or centWave and you are sure your input data are
# centroided, can ignore warning message "It looks like this file is in profile
# mode. [method] can process only centroid mode data !" since this is just based
# on a heuristic. That is, you can ignore the message if you are certain data
# are in centroid mode. You can verify this by opening one of your converted
# .mzXML files in a text reader. You should see:
# <dataProcessing centroided="1"></dataProcessing> (a "0" is bad)
#
# For more on this error, see
# http://metabolomics-forum.com/viewtopic.php?f=8&t=267 or
# https://groups.google.com/forum/#!topic/xcms/xybDDQTaQiY
#
# 2. So long as the number of peak data insertion problems is relatively low
# (i.e., < 100), you can safely ignore the error. Otherwise, might try lowering
# the ppm
#
#  3. On-the-fly plotting features (i.e., with sleep ≥ 0.001 enabled) don't
# appear to function properly in Mac RStudio

#####################################################################################
##### Grouping and retention time correction using xcms (and IPO, if desired) #######
#####################################################################################

if (groupretcor.prams.source==1) { # user wants to run IPO now

    ################# Use IPO to optimize some group, retcor parameters #############

    # define ranges of parameters to be tested

    retcorGroupParameters = getDefaultRetGroupStartingParams(retcorMethod=retcor.meth) # get defaults

    # set some parameter ranges invididually for group.density

    retcorGroupParameters$bw = c(3,15)
    retcorGroupParameters$minfrac = c(0.2,0.5)
    retcorGroupParameters$minsamp = 2
    retcorGroupParameters$mzwid = c(0.001,0.035)
    retcorGroupParameters$profStep = c(0.01,1)

    if (retcor.meth=="loess") {

        # set some parameter ranges invididually for retcor.loess, if retcor.loess was selected
        retcorGroupParameters$missing = c(1,3)
        retcorGroupParameters$extra = c(1,3)
        retcorGroupParameters$smooth = "loess"
        retcorGroupParameters$span = c(0.1,0.3)
        retcorGroupParameters$family = "gaussian" # want to leave outliers in for the time being
        retcorGroupParameters$plottype = "none"
    }

    # perform optimization

    print(paste0("Using R package IPO to optimize group() and retcor() settings with starting parameters user has specified in script."))
    print(paste0("retcor method '",retcor.meth,"' will be used..."))

    resultRetcorGroup = optimizeRetGroup(xset=xset_centWave, params=retcorGroupParameters,
    nSlaves=4, subdir="rsmDirectory")

    ################# Export IPO starting value(s) and optimal settings for each parameter to .csv #############

    # write 3-column table to .csv using write.table()

    if (retcor.meth=="obiwarp") {
        # have to remove resultRetcorGroup$best_settings$center and append it to
        # the end of the concatenated matrix, since there's no option to specify
        # it in retcorGroupParameters

        retcorGroup.exportmat = cbind(sort(rownames(as.matrix(retcorGroupParameters))),
        as.character(resultRetcorGroup$best_settings[-length(resultRetcorGroup$best_settings)][sort(names(resultRetcorGroup$best_settings[-length(resultRetcorGroup$best_settings)]))]),
        as.character(retcorGroupParameters[sort(rownames(as.matrix(retcorGroupParameters)))]))

        retcorGroup.exportmat = rbind(retcorGroup.exportmat,c("center",resultRetcorGroup$best_settings$center,"NA"))

    } else if (retcor.meth=="loess") { # the same right now as in retcor.meth = "obiwarp," but leaving the option here in case we want something different in the future

        # have to remove resultRetcorGroup$best_settings$center and append it to
        # the end of the concatenated matrix, since there's no option to specify
        # it in retcorGroupParameters

        retcorGroup.exportmat = cbind(sort(rownames(as.matrix(retcorGroupParameters))),
        as.character(resultRetcorGroup$best_settings[-length(resultRetcorGroup$best_settings)][sort(names(resultRetcorGroup$best_settings[-length(resultRetcorGroup$best_settings)]))]),
        as.character(retcorGroupParameters[sort(rownames(as.matrix(retcorGroupParameters)))]))

        retcorGroup.exportmat = rbind(retcorGroup.exportmat,c("center",resultRetcorGroup$best_settings$center,"NA"))
    }

    timestamp.now = print(genTimeStamp())

    write.table(retcorGroup.exportmat,
    file = paste("IPO_retcorGroupparamfits_",timestamp.now,".csv",sep=""),
    col.names = c("retcor_or_group_parameter","IPO_optim_value","Starting_value(s)"),
    row.names = FALSE,
    sep=",")

    print(paste0("IPO optimization complete. Optimized and starting values for"))
    print(paste0("group.density and retcor.",retcor.meth," parameters written to file:"))
    print(paste0("IPO_retcorGroupparamfits_",timestamp.now,".csv"))

    # use the just-obtained optimized parameter values for grouping and
    # retention time correction

    print(paste0("Using IPO-optimized settings for group and retcor..."))

    # settings for group.density

    density.bw = resultRetcorGroup$best_settings$bw
    density.max = resultRetcorGroup$best_settings$max
    density.minfrac = resultRetcorGroup$best_settings$minfrac
    density.minsamp = resultRetcorGroup$best_settings$minsamp
    density.mzwid = resultRetcorGroup$best_settings$mzwid

    # settings for selected retcor method

    if (retcor.meth=="obiwarp") {

        obiwarp.profStep = resultRetcorGroup$best_settings$profStep
        obiwarp.response = resultRetcorGroup$best_settings$response
        obiwarp.distFunc = resultRetcorGroup$best_settings$distFunc
        obiwarp.gapInit = resultRetcorGroup$best_settings$gapInit
        obiwarp.gapExtend = resultRetcorGroup$best_settings$gapExtend
        obiwarp.factorDiag = resultRetcorGroup$best_settings$factorDiag
        obiwarp.factorGap = resultRetcorGroup$best_settings$factorGap
        obiwarp.localAlignment = resultRetcorGroup$best_settings$localAlignment

    } else if (retcor.meth=="loess") {

        loess.missing = resultRetcorGroup$best_settings$missing
        loess.extra = resultRetcorGroup$best_settings$extra
        loess.smoothing = resultRetcorGroup$best_settings$smooth
        loess.span = resultRetcorGroup$best_settings$span
        loess.family = resultRetcorGroup$best_settings$family

    }

} else if (groupretcor.prams.source==2) { # user wants to use settings specified below

    print(paste0("Using values of group and retcor parameters specified in the script by user..."))

    # retcor.loess settings below are the function defaults

    loess.missing = 1
    loess.extra = 1
    loess.smoothing = "loess"
    loess.span = c(0.2)
    loess.family = "gaussian" # want to leave outliers in for the time being

    # retcor.obiwarp settings below are the function defaults

    obiwarp.center = NULL
    obiwarp.profStep = 1
    obiwarp.response = 1
    obiwarp.distFunc = "cor_opt"
    obiwarp.gapInit = NULL
    obiwarp.gapExtend = NULL
    obiwarp.factorDiag = 2
    obiwarp.factorGap = 1
    obiwarp.localAlignment = 0

    # settings for group.density below are based on the recommended HPLC/Orbitrap settings from Table 1 of Patti et al., 2012, "Meta-analysis of untargeted metabolomic data from multiple profiling experiment," Nature Protocols 7: 508-516

    density.bw = 5 # 15?
    density.max = 50
    density.minfrac = 0.25
    density.minsamp = 2
    density.mzwid = 0.015 # 0.001?

} else if (groupretcor.prams.source==3) { # user wants to read in parameter values from file

    print(paste0("Loading values of group and retcor parameters from previous IPO optimization run in .csv file:"))
    print(paste0(saved_IPO_params_groupretcor))

    groupretcorprams.from.file = read.csv(saved_IPO_params_groupretcor,colClasses = "character")

    # read in group.density parameter values

    density.bw = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="bw",2])
    density.max = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="max",2])
    density.minfrac = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="minfrac",2])
    density.minsamp = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="minsamp",2])
    density.mzwid = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="mzwid",2])

    # read in parameter values for selected retcor method

    if (retcor.meth=="obiwarp") {

        obiwarp.profStep = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="profStep",2])
        obiwarp.response = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="response",2])
        obiwarp.distFunc = groupretcorprams.from.file[groupretcorprams.from.file[,1]=="distFunc",2]
        obiwarp.gapInit = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="gapInit",2])
        obiwarp.gapExtend = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="gapExtend",2])
        obiwarp.factorDiag = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="factorDiag",2])
        obiwarp.factorGap = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="factorGap",2])
        obiwarp.localAlignment = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="localAlignment",2])

    } else if (retcor.meth=="loess") {

        loess.missing = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="missing",2])
        loess.extra = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="extra",2])
        loess.smoothing = groupretcorprams.from.file[groupretcorprams.from.file[,1]=="smooth",2]
        loess.span = as.numeric(groupretcorprams.from.file[groupretcorprams.from.file[,1]=="span",2])
        loess.family = groupretcorprams.from.file[groupretcorprams.from.file[,1]=="family",2]

    }

}

# specify some additional settings we wish to keep constant, regardless of where the parameters above were obtained

obiwarp.center = NULL
obiwarp.plottype = "deviation" # "none"
density.sleep = 0
loess.plottype = "mdevden" # none

################# Perform grouping and retention time correction on dataset #############

print(paste0("Performing grouping and retention time correction on dataset"))
print(paste0("Using group.density and retcor.",retcor.meth))

# initial grouping

# # method "nearest" with settings below seems to work better than method = "density," but takes absolutely forever; however, it seems to take less time crunching centWave picked data than massifquant picked data

# xset_centWave = group(xset_centWave,
# method = "nearest",
# mzVsRTbalance=10,
# mzCheck=0.2,
# rtCheck=30,
# kNN=10
# )

# using method = "density" with settings from above

xset_gr = group(xset_centWave,
method = "density",
bw = density.bw,
minfrac = density.minfrac,
minsamp = density.minsamp,
mzwid = density.mzwid,
max = density.max,
sleep = density.sleep)

# chromatographic alignment (retention time correction)

if (retcor.meth=="loess") {
    xset_gr.ret = retcor(xset_gr,
    # method = "loess", # this appears unnecessary
    missing = loess.missing,
    extra = loess.extra,
    smooth = "loess",
    span = loess.span,
    family = loess.family,
    plottype = loess.plottype,
    col = NULL,
    ty = NULL
    )

} else if (retcor.meth=="obiwarp") {

    xset_gr.ret = retcor.peakgroups(xset_gr,
    method = "obiwarp",
    plottype = obiwarp.plottype,
    profStep = obiwarp.profStep,
    center = obiwarp.center,
    response = obiwarp.response,
    distFunc = obiwarp.distFunc,
    gapInit = obiwarp.gapInit,
    gapExtend = obiwarp.gapInit,
    factorDiag = obiwarp.factorDiag,
    factorGap = obiwarp.factorGap,
    localAlignment = obiwarp.localAlignment,
    initPenalty = 0
    )

}

# perform grouping again

print(paste0("Performing second peak grouping after application of retcor..."))

# using method = "density" with settings from above

xset_gr.ret.rg = group(xset_gr.ret,
method = "density",
bw = density.bw,
minfrac = density.minfrac,
minsamp = density.minsamp,
mzwid = density.mzwid,
max = density.max,
sleep = density.sleep
)

# fill missing peaks

print(paste0("Filling missing peaks..."))

xset_gr.ret.rg.fill = fillPeaks.chrom(xset_gr.ret.rg, nSlaves = 4)

#####################################################################################
##### Isotope peak identification, creation of xsAnnotate object using CAMERA #######
#####################################################################################

print(paste0("Applying CAMERA to identify isotopic peaks, create xsAnnotate object, and create CAMERA pseudospectra using correlation of xcms peak groups between and within samples. These pseudospectra are the groups within which the adduct hierarchy and retention time screening criteria will be applied using LOBSTAHS"))

# first, a necessary workaround to avoid a import error; see https://support.bioconductor.org/p/69414/
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

# create annotated xset using wrapper annotate(), allowing us to perform all CAMERA tasks at once

xset_a = annotate(xset_gr.ret.rg.fill,

quick=FALSE, # set to FALSE because we want to run groupCorr; will also cause CAMERA to run adduct annotation. while LOBSTAHS will do its own adduct identification later, it doesn't hurt to do this now if it lets CAMERA create better pseudospectra
sample=NA, # use all samples
nSlaves=4, # use 4 sockets

# group FWHM settings
# using defaults for now

sigma=6,
perfwhm=0.6,

# groupCorr settings
# using defaults for now

cor_eic_th=0.75,
graphMethod="hcs",
pval=0.05,
calcCiS=TRUE,
calcIso=TRUE,
calcCaS=FALSE, # weird results with this set to TRUE

# findIsotopes settings

maxcharge=4,
maxiso=4,
minfrac=0.5, # 0.25?

# adduct annotation settings

psg_list=NULL,
rules=NULL,
polarity=subset.polarity,
multiplier=3,
max_peaks=100,

# common to multiple tasks

intval="into",
ppm=2.5,
mzabs=0.0015

)

cleanParallel(xset_a) # kill sockets

# at this point, should have an xsAnnotate object called "xset_a" in hand, which
# will serve as the primary input to the main screening and annotation function
# "doLOBscreen" in LOBSTAHS

print(paste0("xsAnnotate object 'xset_a' has been created. User can now use LOBSTAHS to perform screening..."))

print(xset_a)
