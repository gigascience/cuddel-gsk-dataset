library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)

#### Data import ####

## Get the full path to the CDF files
cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)

## Create a phenodata data.frame to describe experimental setup
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF", replacement = "", fixed = TRUE),
    sample_group = c(rep("KO", 6), rep("WT", 6)),
    stringsAsFactors = FALSE)

# load the raw data as an OnDiskMSnExp object
raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

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
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")

## Plot all chromatograms - this is the total ion chromatogram
pdf("tic_plot.pdf")
plot(bpis, col = group_colors[raw_data$sample_group])
dev.off()

# Extract chromatogram of first sample, check retention times
bpi_1 <- bpis[1, 1]
head(rtime(bpi_1))

# Check intensity values
head(intensity(bpi_1))

# Create boxplots representing the distribution of total ion currents per file.
# Such plots can be very useful to spot problematic or failing MS runs.
tc <- split(tic(raw_data), f = fromFile(raw_data))
pdf("boxplot.pdf")
boxplot(tc, col = group_colors[raw_data$sample_group], ylab = "intensity", main = "Total ion current")
dev.off()

#### Chromatographic peak detection ####

# Evaluate a typical chromatographic peak width by plotting EIC for one peak
# Define the rt and m/z range of the peak area
rtr <- c(2700, 2900)
mzr <- c(334.9, 335.1)
# Extract chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
pdf("peak.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

# Extract full MS data (intensity, retention time and m/z values) of above peak
raw_data %>%
    filterRt(rt = rtr) %>%
    filterMz(mz = mzr) %>%
    plot(type = "XIC")

# Detect peak
cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
xdata <- findChromPeaks(raw_data, param = cwp, BPPARAM = bpparam())

# Look at peak detection results
head(chromPeaks(xdata))

# Calculate per-file summaries using above data
summary_fun <- function(z) {
    c(peak_count = nrow(z),
    rt = quantile(z[, "rtmax"] - z[, "rtmin"]))}
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
# space for 3rd sample
pdf("sample3peaks.pdf")
plotChromPeaks(xdata, file = 3)
dev.off()


# Plot frequency of identified peaks per file along the retention time axis.
# This allows to identify time periods along the MS run in which a higher number
# of peaks was identified and evaluate whether this is consistent across files
pdf("peak_freq.pdf")
plotChromPeakImage(xdata)
dev.off()

# Highlight identified chromatographic peaks for the example peak from before.
# Evaluating such plots on a list of peaks corresponding to known peaks or
# internal standards helps to ensure that peak detection settings were
# appropriate and correctly identified the expected peaks.
pdf("highlight_peak.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata,
    border = group_colors[chr_raw$sample_group],
    lty = 3,
    rt = rtr,
    mz = mzr,
    type = "rect")
dev.off()

# Draw each peak as a polygon
pdf("peak_polygons.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata, col = group_colors[chr_raw$sample_group],
    lty = 3, rt = rtr, mz = mzr, border = NA, type = "polygon")
dev.off()

# extract identified chromatographic peaks for a selected region by providing
# the respective m/z and retention time ranges with the mz and rt arguments in
# the chromPeaks method
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
    caption = paste("Identified chromatographic peaks in a selected ",
    "m/z and retention time range."))

# Plot also the distribution of peak intensity per sample. This allows to
# investigate whether systematic differences in peak signals between samples are
# present.
## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]), f = chromPeaks(xdata)[, "sample"])
pdf("peak_intensity_distribution.pdf")
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
    ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL)
dev.off()

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
pdf("BPC_adjustedRT.pdf")
plot(bpis_adj, col = group_colors[bpis_adj$sample_group])
dev.off()

## Plot also the difference of adjusted to raw retention time.
# Too large differences between adjusted and raw retention times could indicate
# poorly performing samples or alignment.
pdf("BPC_diff_adjusted_rawRT.pdf")
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
pdf("BPC_diff_adjusted_rawRT.pdf")
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group],
    peakGroupsCol = "grey", peakGroupsPch = 1)
dev.off()

# At last we evaluate the impact of the alignment on the test peak.
pdf("single_peak.pdf")
par(mfrow = c(2, 1))
## Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
pdf("single_peak_RTadjusted.pdf")
plot(chr_adj, col = group_colors[chr_raw$sample_group])
dev.off()


########################
#### Correspondence ####
########################

# Correspondence matches peaks between samples using groupChromPeaks function.
# The peak density method [5] is used to group chromatographic peaks by
# combining chromatographic peaks depending on the density of peaks along the
# retention time axis within small slices along the mz dimension

# plot chromatogram for an mz slice with multiple chromatographic peaks within
# each sample. We use below a value of 0.4 for the minFraction parameter hence
# only chromatographic peaks present in at least 40% of the samples per sample
# group are grouped into a feature. The sample group assignment is specified
# with the sampleGroups argument.
## Define the mz slice.
mzr <- c(305.05, 305.15)

## Extract and plot the chromatograms
chr_mzr <- chromatogram(xdata, mz = mzr, rt = c(2500, 4000))
par(mfrow = c(3, 1), mar = c(1, 4, 1, 0.5))
cols <- group_colors[chr_mzr$sample_group]
plot(chr_mzr, col = cols, xaxt = "n", xlab = "")
## Highlight the detected peaks in that region.
highlightChromPeaks(xdata, mz = mzr, col = cols, type = "point", pch = 16)
## Define the parameters for the peak density method
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
    minFraction = 0.4, bw = 30)
par(mar = c(4, 4, 1, 0.5))
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
    pch = 16, xlim = c(2500, 4000))
## Use a different bw
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
    minFraction = 0.4, bw = 20)
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
    pch = 16, xlim = c(2500, 4000))

