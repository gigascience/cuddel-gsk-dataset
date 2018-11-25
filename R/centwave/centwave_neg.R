# Title     : centwave_neg.R
# Objective : TODO
# Created by: peterli
# Created on: 16/11/2018

# https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html

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

# Look at retention time values
head(rtime(raw_data))

# Organise mz values by file
mzs <- mz(raw_data)
mzs_by_file <- split(mzs, f = fromFile(raw_data))  # Split the list by file
length(mzs_by_file)

# Plot base peak chromatogram for each file. aggregationFun is set to "max" to
# return for each spectrum the maximal intensity thereby creating the BPC from
# raw data
bpis <- chromatogram(raw_data, aggregationFun = "max")
# Define colors for the 4 block groups
group_colors <- paste0(brewer.pal(4, "Set1")[1:4], "60")
names(group_colors) <- c("block1", "block2", "block3", "block4")

## Plot all chromatograms from the files
pdf("./centwave/bpc_plot.pdf", width = 20, height = 8)
plot(bpis, col = group_colors[raw_data$sample_group])
dev.off()

## Get the total ion chromatogram. This reads data from the files.
tis <- chromatogram(raw_data, aggregationFun = "sum")
## Define colors for the two groups
group_colors <- paste0(brewer.pal(4, "Set1")[1:4], "60")
names(group_colors) <- c("block1", "block2", "block3", "block4")

## Plot all chromatograms - this is the total ion chromatogram
pdf("./centwave/tic_plot.pdf", width = 20, height = 8)
plot(tis, col = group_colors[raw_data$sample_group])
dev.off()

# Extract chromatogram of first sample
bpi_1 <- bpis[1, 1]
# Access its retention time
head(rtime(bpi_1))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##   0.3672   0.5241   1.1115   1.7150   2.3154   2.8934
# Access its intensity values
head(intensity(bpi_1))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##     8278  8567647  7064262  6896613  9335845  6393502

# Create boxplots representing distribution of total ion currents per file.
# These plots are useful to spot problematic or failing MS runs.
tc <- split(tic(raw_data), f = fromFile(raw_data))
pdf("./centwave/boxplot.pdf")
boxplot(tc, col = group_colors[raw_data$sample_group], ylab = "intensity", main = "Total ion current")
dev.off()


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