# Title     : centwave_neg.R
# Objective : Explore mass spectra using xcms and other mass spec R packages
# Created by: peterli
# Created on: 16/11/2018

# https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html

library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)

# cd /Users/peterli/PhpstormProjects/cuddel-gsk-dataset/R/centwave

# Helper functions for selecting data files
source('../functions.R')

#### Data import ####

# Location of GSK data set
datadir = "/Users/peterli"
# File path for negative files
neg_dir = paste(datadir, "/gsk/cdf/esi_neg/netcdf", sep="")
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

# Get 2 QC samples from each block to test centwave processing
test_qcs <- getQCSampleNamesByBlock(block="1")[1:2]
test_qcs <- append(test_qcs, getQCSampleNamesByBlock(block="2")[1:2])
test_qcs <- append(test_qcs, getQCSampleNamesByBlock(block="3")[1:2])
test_qcs <- append(test_qcs, getQCSampleNamesByBlock(block="4")[1:2])
test_qcs_paths <- lapply(test_qcs, getGSKFilePath)
test_qcs_paths <- unlist(test_qcs_paths)

# Create a phenodata data.frame to describe batch setup
pd <- data.frame(sample_name = sub(basename(test_qcs_paths), pattern = ".cdf", replacement = "", fixed = TRUE),
    sample_group = c(rep("block1", 2), rep("block2", 2), rep("block3", 2), rep("block4", 2)),
    stringsAsFactors = FALSE)

# Load raw data as OnDiskMSnExp object
raw_data <- readMSData(files = test_qcs_paths, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

#### Initial data inspection ####

# Get total number of spectra
length(raw_data)
## [1] 18423

# Look at retention time values
head(rtime(raw_data))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##   0.3672   0.5241   1.1115   1.7150   2.3154   2.8934

# Number of RTs
length(rtime(raw_data))
## [1] 18423
# Basically there are spectra from 18,423 retention times

# Organise mz values by file
mzs <- mz(raw_data)
mzs_by_file <- split(mzs, f = fromFile(raw_data))  # Split the list by file
length(mzs_by_file)
## [1] 8

# LC/MS can be used to identify and quantify metabolites. They are separated
# within a column and, at each time point when they come out of the column, the
# abundance (intensity) of a particular group of metabolites is measured by the
# mass spectrometer. This is called a chromatogram.

# MS data can be read from the files to get the base peak chromatogram. This
# chromatogram monitors only the most intense peak in each spectrum. This means
# that the base peak chromatogram represents the intensity of the most intense
# peak at every point in the analysis. Base peak chromatograms often have a
# cleaner look and thus are more informative than TIC chromatograms because the
# background is reduced by focusing on a single analyte at every point.
bpis <- chromatogram(raw_data, aggregationFun = "max")
# Define colors for the 4 block groups
group_colors <- paste0(brewer.pal(4, "Set1")[1:4], "60")
names(group_colors) <- c("block1", "block2", "block3", "block4")

## Plot all chromatograms from the files
pdf("./output/gsk/bpc_plot.pdf", width = 20, height = 8)
plot(bpis, col = group_colors[raw_data$sample_group])
dev.off()

# Plot all chromatograms - this is the total ion chromatogram or TIC and
# represents the summed intensity across the entire range of masses being
# detected at every point in the analysis. The range is typically several
# hundred mass-to-charge units or more. In complex samples, the TIC chromatogram
# often provides limited information as multiple analytes elute simultaneously,
# obscuring individual species.
tis <- chromatogram(raw_data, aggregationFun = "sum")
## Define colors for the 4 block groups
group_colors <- paste0(brewer.pal(4, "Set1")[1:4], "60")
names(group_colors) <- c("block1", "block2", "block3", "block4")

## Plot all chromatograms - this is the total ion chromatogram
pdf("./output/gsk/tic_plot.pdf", width = 20, height = 8)
plot(tis, col = group_colors[raw_data$sample_group])
dev.off()

# Extract chromatogram of first sample, check retention times
bpi_1 <- bpis[1, 1]
# Access its retention time
head(rtime(bpi_1))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##   0.3672   0.5241   1.1115   1.7150   2.3154   2.8934
# Check intensity values
head(intensity(bpi_1))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##     8278  8567647  7064262  6896613  9335845  6393502

# Create boxplots representing the distribution of total ion currents per file.
# Such plots can be very useful to spot problematic or failing MS runs.
tc <- split(tic(raw_data), f = fromFile(raw_data))
pdf("./output/gsk/boxplot.pdf")
boxplot(tc, main="Distribution of total ion currents from GSK dataset files",
col = group_colors[raw_data$sample_group],
ylab = "intensity", main = "Total ion current")
dev.off()

# We can cluster the samples based on similarity of their base peak
# chromatogram. This can also be helpful to spot potentially problematic samples
# in an experiment or generally get an initial overview of the sample grouping
# in the experiment. Since retention times between samples are not exactly
# identical, the bin function is used to group intensities in fixed time ranges
# (bins) along the retention time axis. In the present example, we use a bin
# size of 1 second, the default is 0.5 seconds. Clustering is performed using
# complete linkage hierarchical clustering on the pairwise correlations of the
# binned base peak chromatograms.

## Bin the BPC
bpis_bin <- bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name

## Perform the cluster analysis
pdf("./output/gsk/bpc_similarity.pdf")
pheatmap(cormat, annotation = ann, annotation_color = list(group = group_colors))
dev.off()

#### Chromatographic peak detection ####

# At each retention time point, a group of metabolites is ionized (charged) and
# fired through a mass spectrometer. The metabolite ions can be separated based
# on their mass-to-charge ratio. This adds another dimension to the graph: an
# axis with the mass-to-charge ratio of the metabolites found at that time
# point, and their intensities.

# Therefore, there are multiple spectra. Each spectrum can be simplified into a
# graph of peaks (local maxima).

# To get m/z against intensity plot for a specific RT

# List spectra
spectra(raw_data)

# Get the number of spectra present
length(spectra(raw_data))
## [1] 18423

# Get second spectrum.
sp2 <- raw_data[[2]]

# Get RT for second spectrum
rtime(sp2)
## [1] 0.5241

# Plot second spectrum which is a graph of m/z against signal intensity
pdf("./output/gsk/sp2.pdf")
plot(sp2)
dev.off()

# Define mz and RT ranges to extract main peak from second spectrum
pdf("./output/gsk/peak.pdf")
peak <- filterMz(filterRt(raw_data, rt = c(0.5241, 0.5241)), mz = c(50, 125))
# XIC displays a combined plot of RT vs m/z and RT vs largest signal in spectrum
plot(peak, xaxt="none", type="XIC")
axis(1, seq(0.3, 0.7, 0.05))
dev.off()

# Zoom into a part of the above spectrum using above m/z and RT ranges
pdf("./output/gsk/zoom_peak.pdf")
plot(peak, centroided = TRUE)
dev.off()

#### Find tryptophan peak ####

# To look for the peak corresponding to tryptophan, we need to get the molecular
# mass of the amino acid and subtract mass of 1 electron from it. Then look for
# the peak at this m/z value in negative ionisation mode data.
#
# Used this [mass calculator](https://www.envipat.eawag.ch/index.php)

# Using the molecular formula of Tryptophan C11H12N2O2, its molecular mass is
# 204.0899 mol mass without adducts - use 3 decimal points for orbitrap ppm.

# Tryptophan as a negative ion has a monoisotopic mass 203.0826012.

# For tryptophan as a positive ion from the extra proton, the monoisotopic mass
# is 205.0971541. Use this m/z value to look for the tryptophan peak in positive
# ion mode data

ht(rtime(raw_data))
F1.S0001  F1.S0002  F1.S0003 F1.S0004  F1.S0005  F1.S0006
[1,]    0.3672    0.5241    1.1115    1.715    2.3154    2.8934
[2,] 1436.5604 1437.1133 1437.7089 1438.290 1438.8683 1439.4577


## Define the rt and m/z range of the peak area
rtr <- c(0.3672, 1439.4577)
mzr <- c(203.0820, 203.083)
## extract the chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
pdf("./output/gsk/peak_tryptophan.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

peak <- filterMz(filterRt(raw_data, rt = c(0.5241, 0.5241)), mz = c(50, 125))
pdf("./output/gsk/peak_tryptophan.pdf")
plot(peak, xaxt="none", type="XIC")
dev.off()


#### Handling raw MS data ####

# The mzR package provides an interface to the proteowizard C/C++ code base to
# access various raw data files, such as mzML, mzXML, netCDF, and mzData. The
# data is accessed on-disk, i.e it is not loaded entirely in memory by default
# but only when explicitly requested. The three main functions are openMSfile to
# create a file handle to a raw data file, header to extract metadata about the
# spectra contained in the file and peaks to extract one or multiple spectra of
# interest. Other functions such as instrumentInfo, or runInfo can be used to
# gather general information about a run.
#
# Below, we access the raw data file downloaded in the previous section and open
# a file handle that will allow us to extract data and metadata of interest.

# Get metadata about spectra contained in files using header function
hd <- header(raw_data)

# Get number of scans or spectra. There are 18423 scans/spectra
dim(hd)
## [1] 18423    31

names(hd)
## [1] "fileIdx"                  "spIdx"
## [3] "smoothed"                 "seqNum"
## [5] "acquisitionNum"           "msLevel"
## [7] "originalPeaksCount"       "totIonCurrent"
## [9] "retentionTime"            "basePeakMZ"
## [11] "basePeakIntensity"        "collisionEnergy"
## [13] "ionisationEnergy"         "highMZ"
## [15] "precursorScanNum"         "precursorMZ"
## [17] "precursorCharge"          "precursorIntensity"
## [19] "mergedScan"               "mergedResultScanNum"
## [21] "mergedResultStartScanNum" "mergedResultEndScanNum"
## [23] "injectionTime"            "spectrumId"
## [25] "centroided"               "ionMobilityDriftTime"
## [27] "polarity"                 "spectrum"


# Extract metadata and scan data for scan 1000
hd[1000, ]
##          fileIdx spIdx smoothed seqNum acquisitionNum msLevel
## F1.S1000       1  1000       NA   1000           1000       1
##          originalPeaksCount totIonCurrent retentionTime basePeakMZ
## F1.S1000                  1      20235066         557.4         -1
##          basePeakIntensity collisionEnergy ionisationEnergy highMZ
## F1.S1000                -1              -1               -1     -1
##          precursorScanNum precursorMZ precursorCharge precursorIntensity
## F1.S1000               -1          -1              -1                 -1
##          mergedScan mergedResultScanNum mergedResultStartScanNum
## F1.S1000         -1                  -1                       -1
##          mergedResultEndScanNum injectionTime spectrumId centroided
## F1.S1000                     -1            -1  scan=1000         NA
##          ionMobilityDriftTime polarity spectrum
## F1.S1000                   -1       NA     1000

# Extract spectra from scan 1000
head(peaks(raw_data, 1000))
## [,1]     [,2]
## [1,] 52.33684 3490.230
## [2,] 52.63135 3543.195
## [3,] 53.25687 3011.438q
## [4,] 54.01475 3457.050
## [5,] 54.73533 3067.098
## [6,] 55.42279 3363.786

# Plot peaks from scan 1000
plot(peaks(ms, 1000), type = "h")

# Plot specific set of mass spectra using mzR functions
ms1 <- which(hd$msLevel == 1)
# Specify we want spectra from compounds eluted between 30-35 minutes RT
rtsel <- hd$retentionTime[ms1] / 60 > 20 & hd$retentionTime[ms1] / 60 < 22

# The map
M <- MSmap(ms, ms1[rtsel], 52, 53, .005, hd, zeroIsNA = TRUE)

# Plot slice of data
plot(M, aspect = 1, allTicks = FALSE)

# Plot 3D graph of data slice
plot3D(M)

########################################
#### Chromatographic peak detection ####
########################################

# Look at a single peak for a compound to check centWave parameter settings.
# The peak looks strange with our orbitrap data.
rtr <- c(34.0, 46.0)
mzr <- c(238.5, 238.9)
# Extract the chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
pdf("./centwave/peak.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

raw_data %>%
    filterRt(rt = rtr) %>%
    filterMz(mz = mzr) %>%
    plot(type = "XIC")

# Perform chromatographic peak detection using centWave algorithm. Noise
# parameter is set to 5000 to reduce run time
cwp <- CentWaveParam(ppm=20, peakwidth = c(7, 25), noise = 5000)
xdata <- findChromPeaks(raw_data, param = cwp, BPPARAM = bpparam())

# Access results with chromPeaks method. The returned matrix provides the m/z
# and retention time range for each identified chromatographic peak as well as
# the integrated signal intensity (“into”) and the maximal peak intensitity
# (“maxo”). Columns “sample” contains the index of the sample in the
# object/experiment in which the peak was identified.
head(chromPeaks(xdata))

# Summarise results per file
summary_fun <- function(z) {
    c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))}
T <- lapply(split.data.frame(chromPeaks(xdata),
    f = chromPeaks(xdata)[, "sample"]),
    FUN = summary_fun)
T <- do.call(rbind, T)
rownames(T) <- basename(fileNames(xdata))
pandoc.table(T,
    caption = paste0("Summary statistics on identified chromatographic",
    " peaks. Shown are number of identified peaks per",
    " sample and widths/duration of chromatographic ",
    "peaks."))

# Plot location of identified chromatographic peaks in the m/z - retention time
# space for one file using the plotChromPeaks function for 3rd sample
pdf("./centwave/location_peaks_3rdsample.pdf")
plotChromPeaks(xdata, file = 3)
dev.off()

# Get global overview of the peak detection by plotting frequency of identified
# peaks per file along the retention time axis. This allows to
# identify time periods along the MS run in which a higher number of peaks was
# identified and evaluate whether this is consistent across files.
pdf("./centwave/peak_freq_per_file.pdf")
plotChromPeakImage(xdata)
dev.off()

# Highlight identified chromatographic peaks for example peak from before.
# Evaluating such plots on a list of peaks corresponding to known peaks or
# internal standards helps to ensure that peak detection settings were
# appropriate and correctly identified the expected peaks.
pdf("./centwave/highlight_peak.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata,
    border = group_colors[chr_raw$sample_group],
    lty = 3, rt = rtr, mz = mzr, type = "rect")
dev.off()

# specifically extract identified chromatographic peaks for a selected region by
# providing the respective m/z and retention time ranges with the mz and rt
# arguments in the chromPeaks method.
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
caption = paste("Identified chromatographic peaks in a selected ",
    "m/z and retention time range."))

# Plot distribution of peak intensity per sample to check whether systematic
# differences in peak signals between samples are present.
pdf("./centwave/peak_intensity_distribution.pdf")
# Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]), f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
    ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL)
dev.off()


###################
#### Alignment ####
###################

# Analytes can elute at different times between samples during chromatography.
# Retention time correction can align same peak signals between different
# samples within an experiment. Use adjustRtime function to perform RT
# correction to align signals between samples.
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

# Extract adjusted retention times
head(adjustedRtime(xdata))

# Plot BPC on adjusted data to evaluate impact of signal alignment using RT
# Result looks strange and some runs look truncated
bpis_adj <- chromatogram(xdata, aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
pdf("./centwave/BPC_adjustedRT.pdf", width = 20, height = 8)
plot(bpis_adj, col = group_colors[bpis_adj$sample_group])
dev.off()

# Plot difference of adjusted to raw retention time. Large differences between
# adjusted and raw retention times could indicate poorly performing samples or
# alignment
pdf("./centwave/BPC_diff_adjusted_rawRT.pdf")
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])
dev.off()

## Does the object have adjusted retention times?
hasAdjustedRtime(xdata)
## Drop the alignment results.
xdata <- dropAdjustedRtime(xdata)

## Does the object have adjusted retention times?
hasAdjustedRtime(xdata)

## Correspondence: group peaks across samples.
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group, minFraction = 0.8)
xdata <- groupChromPeaks(xdata, param = pdp)

## Now the retention time correction.
pgp <- PeakGroupsParam(minFraction = 0.85)

## Get the peak groups that would be used for alignment.
xdata <- adjustRtime(xdata, param = pgp)

## Plot the difference of adjusted to raw retention time.
pdf("./centwave/BPC_diff_adjusted_rawRT.pdf")
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group],
peakGroupsCol = "grey", peakGroupsPch = 1)
dev.off()

# Evaluate impact of RT alignment on test peak
pdf("./centwave/single_peak.pdf")
par(mfrow = c(2, 1))
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
pdf("./centwave/single_peak_RTadjusted.pdf")
plot(chr_adj, col = group_colors[chr_raw$sample_group])
dev.off()