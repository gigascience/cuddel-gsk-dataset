library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)

# Commands from https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html

#### Data import ####

# Get the full path to the CDF files
cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)

## Create a phenodata data.frame to describe experimental setup
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF", replacement = "", fixed = TRUE),
    sample_group = c(rep("KO", 6), rep("WT", 6)),
    stringsAsFactors = FALSE)

# Load the raw data as an OnDiskMSnExp object
raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

#### Initial data inspection ####

# Look at retention time values
head(rtime(raw_data))
## F01.S0001 F01.S0002 F01.S0003 F01.S0004 F01.S0005 F01.S0006
##  2501.378  2502.943  2504.508  2506.073  2507.638  2509.203

# Organise mz values by file
mzs <- mz(raw_data)
mzs_by_file <- split(mzs, f = fromFile(raw_data))  # Split the list by file
length(mzs_by_file)
## [1] 12

# Get the base peak chromatograms. This reads data from the files. This
# chromatogram monitors only the most intense peak in each spectrum. This means
# that the base peak chromatogram represents the intensity of the most intense
# peak at every point in the analysis. Base peak chromatograms often have a
# cleaner look and thus are more informative than TIC chromatograms because the
# background is reduced by focusing on a single analyte at every point.
bpis <- chromatogram(raw_data, aggregationFun = "max")
# Define colors for the two groups
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")

par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
pdf("output/test/BPC.pdf")
plot(bpis, main = "Base peak chromatogram", col = group_colors[bpis$sample_group])
dev.off()

# Plot all chromatograms - this is the total ion chromatogram or TIC and
# represents the summed intensity across the entire range of masses being
# detected at every point in the analysis. The range is typically several
# hundred mass-to-charge units or more. In complex samples, the TIC chromatogram
# often provides limited information as multiple analytes elute simultaneously,
# obscuring individual species.
pdf("output/test/tic_plot.pdf")
plot(bpis, main="Total ion chromatogram from faahKO dataset",
    col = group_colors[raw_data$sample_group])
dev.off()

# Extract chromatogram of first sample, check retention times
bpi_1 <- bpis[1, 1]
head(rtime(bpi_1))
## F01.S0001 F01.S0002 F01.S0003 F01.S0004 F01.S0005 F01.S0006
##  2501.378  2502.943  2504.508  2506.073  2507.638  2509.203

# Check intensity values
head(intensity(bpi_1))
## F01.S0001 F01.S0002 F01.S0003 F01.S0004 F01.S0005 F01.S0006
##     43888     43960     43392     42632     42200     42288

# Create boxplots representing the distribution of total ion currents per file.
# Such plots can be very useful to spot problematic or failing MS runs.
tc <- split(tic(raw_data), f = fromFile(raw_data))
pdf("output/test/boxplot.pdf")
boxplot(tc, main="Distribution of total ion currents from faahKO dataset files",
    col = group_colors[raw_data$sample_group],
    ylab = "intensity", main = "Total ion current")
dev.off()

#### Chromatographic peak detection ####

# Evaluate a typical chromatographic peak width by plotting the extracted ion
# chromatogram (EIC). An EIC displays one or more peaks which represents one or
# more m/z values representing one or more analytes of interest are recovered
# ('extracted') from the entire data set for a chromatographic run.
# Define rt and m/z range of the peak area
rtr <- c(2700, 2900)
mzr <- c(334.9, 335.1)
# Extract chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
pdf("output/test/peak.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

# Extract full MS data (intensity, retention time and m/z values) of above peak
pdf("output/test/XICs.pdf")
raw_data %>%
    filterRt(rt = rtr) %>%
    filterMz(mz = mzr) %>%
    plot(type = "XIC")
dev.off()

# Detect peaks using centwave. Noise is set to 500 reduce running time.
cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000)
xdata <- findChromPeaks(raw_data, param = cwp, BPPARAM = bpparam())

# Results returned in XCMSnExp object and are accessed with chromPeaks method
head(chromPeaks(xdata))
##           mz mzmin mzmax       rt    rtmin    rtmax      into      intb  maxo
## CP0001 453.2 453.2 453.2 2509.203 2501.378 2527.982 1007409.0 1007380.8 38152
## CP0002 236.1 236.1 236.1 2518.593 2501.378 2537.372  253501.0  226896.3 12957
## CP0003 594.0 594.0 594.0 2601.535 2581.191 2637.529  161042.2  149297.3  7850
## CP0004 577.0 577.0 577.0 2604.665 2581.191 2626.574  136105.2  129195.5  6215
## CP0005 369.2 369.2 369.2 2587.451 2556.151 2631.269  483852.3  483777.1  7215
## CP0006 369.2 369.2 369.2 2568.671 2557.716 2578.061  144624.8  144602.9  7033
##           sn sample is_filled
## CP0001 38151      1         0
## CP0002    11      1         0
## CP0003    13      1         0
## CP0004    13      1         0
## CP0005  7214      1         0
## CP0006  7032      1         0
# The returned matrix provides the m/z and retention time range for each
# identified chromatographic peak as well as the integrated signal intensity
# (“into”) and the maximal peak intensitity (“maxo”). Columns “sample” contains
# the index of the sample in the object/experiment in which the peak was
# identified.

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
    " peaks. Shown are number of identified peaks per sample and",
    " widths/duration of chromatographic peaks."))

## ------------------------------------------------------------------------
##     &nbsp;      peak_count   rt.0%   rt.25%   rt.50%   rt.75%   rt.100%
## -------------- ------------ ------- -------- -------- -------- ---------
## **ko15.CDF**      561       7.824   28.17    39.12    50.08     139.3
##
## **ko16.CDF**      741       4.695   23.47    39.12    50.08     148.7
##
## **ko18.CDF**      448       4.694   28.17    42.25    53.21     173.7
##
## **ko19.CDF**      333       7.824    31.3    46.95    59.47     181.5
##
## **ko21.CDF**      276       1.565   35.99    48.51    61.03     164.3
##
## **ko22.CDF**      325       7.824   26.61    43.82    56.34     126.8
##
## **wt15.CDF**      564       4.695    26.6    39.12    50.08     161.2
##
## **wt16.CDF**      507       7.824   28.17    40.69    53.21     151.8
##
## **wt18.CDF**      443       3.13    26.61    42.25    57.12     192.5
##
## **wt19.CDF**      350       7.824    31.3    45.38     57.9     184.7
##
## **wt21.CDF**      317       3.13    28.17    45.38    61.03     172.1
##
## **wt22.CDF**      437       1.565   26.61    40.69    53.21     161.2
## ------------------------------------------------------------------------
##
## Table: Summary statistics on identified chromatographic peaks. Shown are
## number of identified peaks per sample and widths/duration of chromatographic
## peaks.

# Plot location of identified chromatographic peaks in the m/z - retention time
# space for third sample
pdf("output/test/sample3peaks.pdf")
plotChromPeaks(xdata, file = 3)
dev.off()

# Get a global overview of the peak detection by plotting frequency of
# identified peaks per file along the retention time axis.This allows to
# identify time periods along the MS run in which a higher number of peaks was
# identified and evaluate whether this is consistent across files.
pdf("output/test/peak_freq.pdf")
plotChromPeakImage(xdata)
dev.off()

# Highlight identified chromatographic peaks for the example peak from before.
# Evaluating such plots on a list of peaks corresponding to known peaks or
# internal standards helps to ensure that peak detection settings were
# appropriate and correctly identified the expected peaks.
pdf("output/test/highlight_peak.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata,
    border = group_colors[chr_raw$sample_group],
    lty = 3,
    rt = rtr,
    mz = mzr,
    type = "rect")
dev.off()

# Draw each peak as a polygon
pdf("output/test/peak_polygons.pdf")
plot(chr_raw, col = group_colors[chr_raw$sample_group], lwd = 2)
highlightChromPeaks(xdata, col = group_colors[chr_raw$sample_group],
    lty = 3, rt = rtr, mz = mzr, border = NA, type = "polygon")
dev.off()

# Extract identified chromatographic peaks for a selected region by providing
# the respective m/z and retention time ranges with the mz and rt arguments in
# the chromPeaks method
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
    caption = paste("Identified chromatographic peaks in a selected ",
    "m/z and retention time range."))

## -----------------------------------------------------------------------------
##    &nbsp;     mz    mzmin   mzmax    rt    rtmin   rtmax    into      intb
## ------------ ----- ------- ------- ------ ------- ------- --------- ---------
##  **CP0057**   335    335     335    2782   2761    2810    412134    383167
##
##  **CP0628**   335    335     335    2786   2764    2813    1496244   1461187
##
##  **CP1349**   335    335     335    2786   2764    2813    211297    194283
##
##  **CP1795**   335    335     335    2783   2763    2813    285228    264249
##
##  **CP2120**   335    335     335    2797   2775    2816    159059    149230
##
##  **CP3308**   335    335     335    2786   2764    2813    932645    915334
##
##  **CP3838**   335    335     335    2788   2766    2821    1636669   1603951
##
##  **CP4254**   335    335     335    2785   2760    2810    643673    631119
##
##  **CP4601**   335    335     335    2792   2769    2824    876586    848569
##
##  **CP4906**   335    335     335    2789   2774    2807     89583     83927
## -----------------------------------------------------------------------------
##
## Table: Identified chromatographic peaks in a selected  m/z and retention time range. (continued below)
##
##
## ----------------------------------------------
##    &nbsp;     maxo    sn   sample   is_filled
## ------------ ------- ---- -------- -----------
##  **CP0057**   16856   23     1          0
##
##  **CP0628**   58736   72     2          0
##
##  **CP1349**   8158    15     3          0
##
##  **CP1795**   9154    18     4          0
##
##  **CP2120**   6295    13     5          0
##
##  **CP3308**   35856   66     8          0
##
##  **CP3838**   54640   78     9          0
##
##  **CP4254**   20672   54     10         0
##
##  **CP4601**   27200   36     11         0
##
##  **CP4906**   5427    12     12         0
## ----------------------------------------------

# Plot the distribution of peak intensity per sample to investigate whether
# systematic differences in peak signals between samples are present
# Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]), f = chromPeaks(xdata)[, "sample"])
pdf("output/test/peak_intensity_distribution.pdf")
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
    ylab = expression(log[2]~intensity), main = "Peak intensities per sample")
grid(nx = NA, ny = NULL)
dev.off()

###################
#### Alignment ####
###################

# Analytes can elute at different times between samples during chromatography.
# Retention time correction can align same peak signals between different
# samples within an experiment. Use adjustRtime function to perform RT
# correction with obiwarp method
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

## Extract adjusted retention times
head(adjustedRtime(xdata))
## F01.S0001 F01.S0002 F01.S0003 F01.S0004 F01.S0005 F01.S0006
##  2501.378  2502.958  2504.538  2506.118  2507.699  2509.280

# To evaluate the impact of the alignment, plot the BPC on the adjusted data.
# Also plot the differences of the adjusted- to the raw retention times per
# sample using the plotAdjustedRtime function.
# Get the base peak chromatograms.
bpis_adj <- chromatogram(xdata, aggregationFun = "max")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
pdf("output/test/BPC_adjustedRT.pdf")
plot(bpis_adj, main = "Base peak chromatogram after RT correction", col = group_colors[bpis_adj$sample_group])
dev.off()

# Plot the difference of adjusted to raw retention time. Too large differences
# between adjusted and raw retention times could indicate poorly performing
# samples or alignment.
pdf("output/test/BPC_diff_adjusted_rawRT.pdf")
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])
dev.off()

# It is possible to restore the results from the RT correction. Check if the
# object has adjusted retention times
hasAdjustedRtime(xdata)
## [1] TRUE

# Drop RT correction alignment results
xdata <- dropAdjustedRtime(xdata)

# Does the object now have adjusted retention times?
hasAdjustedRtime(xdata)
## [1] FALSE

# We can also use the group peaks method to adjust RT by aligning previously
# identified hook peaks (chromatographic peak features present in most/all
# samples). Ideally, these hook peaks should span most part of the retention
# time range. A first correspondence run is used to identify such peaks

# we can use the adjustedRtimePeakGroups method on the object before alignment
# to evaluate on which features (peak groups) the alignment would be performed.
# This can be useful to test different settings for the peak groups algorithm.
# Also, it is possible to manually select or define certain peak groups (i.e.
# their retention times per sample) and provide this matrix to the
# PeakGroupsParam class with the peakGroupsMatrix argument.

# Default settings are used below, but it is strongly advised to adapt the
# parameters for each data set. The definition of the sample groups (i.e.
# assignment of individual samples to the sample groups in the experiment) is
# mandatory for the PeakDensityParam. If there are no sample groups in the
# experiment, sampleGroups should be set to a single value for each file (e.g.
# rep(1, length(fileNames(xdata))).
# Correspondence: group peaks across samples:
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group, minFraction = 0.8)
xdata <- groupChromPeaks(xdata, param = pdp)

## Now do retention time correction
pgp <- PeakGroupsParam(minFraction = 0.85)

# Get the peak groups that would be used for alignment.
xdata <- adjustRtime(xdata, param = pgp)

# Plot the difference between raw and adjusted retention times using the
# plotAdjustedRtime function, which, if the peak groups method is used for
# alignment, also highlights the peak groups used in the adjustment.
pdf("output/test/BPC_diff_adjusted_rawRT.pdf")
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group],
    peakGroupsCol = "grey", peakGroupsPch = 1)
dev.off()

# We can now evaluate the impact of the alignment on the test peak.
pdf("output/test/single_peak.pdf")
par(mfrow = c(2, 1))
# Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])
dev.off()

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
pdf("output/test/single_peak_RTadjusted.pdf")
plot(chr_adj, col = group_colors[chr_raw$sample_group])
dev.off()


########################
#### Correspondence ####
########################

# Correspondence matches peaks between samples using groupChromPeaks function.
# The peak density method [5] is used to group chromatographic peaks by
# combining chromatographic peaks depending on the density of peaks along the
# retention time axis within small slices along the mz dimension

# Plot a chromatogram for an mz slice with multiple chromatographic peaks within
# each sample. We use below a value of 0.4 for the minFraction parameter hence
# only chromatographic peaks present in at least 40% of the samples per sample
# group are grouped into a feature. The sample group assignment is specified
# with the sampleGroups argument.
# Define the mz slice.
mzr <- c(305.05, 305.15)

## Extract and plot the chromatograms
chr_mzr <- chromatogram(xdata, mz = mzr, rt = c(2500, 4000))
par(mfrow = c(3, 1), mar = c(1, 4, 1, 0.5))
cols <- group_colors[chr_mzr$sample_group]
pdf("output/single_peak_RTadjusted.pdf")
plot(chr_mzr, col = cols, xaxt = "n", xlab = "")
dev.off()
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

