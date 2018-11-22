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

# Extract chromatogram of first sample and access its retention time and
# intensity values
bpi_1 <- bpis[1, 1]
head(rtime(bpi_1))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##   0.3672   0.5241   1.1115   1.7150   2.3154   2.8934

head(intensity(bpi_1))
## F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006
##     8278  8567647  7064262  6896613  9335845  6393502

# Create boxplots representing the distribution of total ion currents per file.
# Such plots can be very useful to spot problematic or failing MS runs
tc <- split(tic(raw_data), f = fromFile(raw_data))
pdf("./centwave/boxplot.pdf")
boxplot(tc, col = group_colors[raw_data$sample_group], ylab = "intensity", main = "Total ion current")
dev.off()


########################################
#### Chromatographic peak detection ####
########################################

# Detect chromatogrpahic peaks with noise set to 500 to reduce run time
cwp <- CentWaveParam(ppm=20, peakwidth = c(7, 25), noise = 5000)
xdata <- findChromPeaks(raw_data, param = cwp, BPPARAM = bpparam())

# results from the chromatographic peak detection can be accessed with the chromPeaks
# method.
head(chromPeaks(xdata))

# Summarise results per file
summary_fun <- function(z) {
    c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
}
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

# plot the location of the identified chromatographic peaks in the m/z -
# retention time space for one file using the plotChromPeaks function. Below we
# plot the chromatographic peaks for the 3rd sample.
plotChromPeaks(xdata, file = 3)

# get a global overview of the peak detection we can plot the frequency of
# identified peaks per file along the retention time axis. This allows to
# identify time periods along the MS run in which a higher number of peaks was
# identified and evaluate whether this is consistent across files.
plotChromPeakImage(xdata)

# highlight the identified chromatographic peaks for the example peak from
# before. Evaluating such plots on a list of peaks corresponding to known peaks
# or internal standards helps to ensure that peak detection settings were
# appropriate and correctly identified the expected peaks.
plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata,
    border = group_colors[chr_raw$sample_group],
    lty = 3, rt = rtr, mz = mzr, type = "rect")

# specifically extract identified chromatographic peaks for a selected region by
# providing the respective m/z and retention time ranges with the mz and rt
# arguments in the chromPeaks method.
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
caption = paste("Identified chromatographic peaks in a selected ",
"m/z and retention time range."))

# plot also the distribution of peak intensity per sample. This allows to
# investigate whether systematic differences in peak signals between samples are
# present.
## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]),
f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL)


###################
#### Alignment ####
###################

# Analytes can elute at different times between samples during chromatography.
# Retention time correction can align same peak signals between different
# samples within an experiment. USe adjestRtime function to perform RT
# correction
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

## Extract adjusted retention times
head(adjustedRtime(xdata))

# To evaluate the impact of the alignment we plot the BPC on the adjusted data
## Get the base peak chromatograms.
bpis_adj <- chromatogram(xdata, aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
pdf("./centwave/BPC_adjustedRT.pdf")
plot(bpis_adj, col = group_colors[bpis_adj$sample_group])
dev.off()

## Plot also the difference of adjusted to raw retention time.
# Too large differences between adjusted and raw retention times could indicate
# poorly performing samples or alignment.
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

# At last we evaluate the impact of the alignment on the test peak.
pdf("./centwave/single_peak.pdf")
par(mfrow = c(2, 1))
## Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
pdf("./centwave/single_peak_RTadjusted.pdf")
plot(chr_adj, col = group_colors[chr_raw$sample_group])
dev.off()