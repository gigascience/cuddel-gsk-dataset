library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)

# XCMS tutorial from https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html

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
# Factors such as changes in pH, temperature (column) and elution solvents
# composition can introduce some variability to chromatographic retention time.
# Retention time correction can correct this shift and align same peak signals
# between different samples within an experiment. Use adjustRtime function to
# perform RT correction with obiwarp method
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

# Extract adjusted retention times
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

# Extract the chromatogram from the adjusted object. You will see that the
# peaks from the different samples are now better aligned after RT adjustment.
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
pdf("output/test/single_peak_RTadjusted.pdf")
plot(chr_adj, col = group_colors[chr_raw$sample_group])
dev.off()


########################
#### Correspondence ####
########################

# Correspondence matches peaks between samples using groupChromPeaks function.
# The peak density method is used to group chromatographic peaks by combining
# chromatographic peaks depending on the density of peaks along the retention
# time axis within small slices along the mz dimension.

# Below, a chromatogram is plotted for an mz slice with multiple chromatographic
# peaks within each sample. A value of 0.4 is used for the minFraction
# parameter so only chromatographic peaks present in at least 40% of the
# samples per sample group are grouped into a feature. The sample group
# assignment is specified with the sampleGroups argument.
# Define the mz slice.
mzr <- c(305.05, 305.15)

# Extract and plot the chromatograms
pdf("output/test/sample_eics.pdf")
chr_mzr <- chromatogram(xdata, mz = mzr, rt = c(2500, 4000))
par(mfrow = c(3, 1), mar = c(1, 4, 1, 0.5))
cols <- group_colors[chr_mzr$sample_group]
plot(chr_mzr, col = cols, xaxt = "n", xlab = "")
dev.off()

# Highlight the detected peaks in that region
highlightChromPeaks(xdata, mz = mzr, col = cols, type = "point", pch = 16)

# Define parameters for peak density method. If bw=30, peaks are combined
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
    minFraction = 0.4, bw = 30)
par(mar = c(4, 4, 1, 0.5))

pdf("output/test/RT_peaks_bw30.pdf")
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
    pch = 16, xlim = c(2500, 4000))
dev.off()

## Using a different lower bw=20 value, the peaks are separated
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
    minFraction = 0.4, bw = 20)
pdf("output/test/RT_peaks_bw20.pdf")
plotChromPeakDensity(xdata, mz = mzr, col = cols, param = pdp,
    pch = 16, xlim = c(2500, 4000))
dev.off()

# Perform correspondence
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
    minFraction = 0.4, bw = 20)
xdata <- groupChromPeaks(xdata, param = pdp)

# Extract feature definition results from correspondence using the
# featureDefinitions method that returns a DataFrame with the definition of the
# features (i.e. the mz and rt ranges and, in column peakidx, the index of the
# chromatographic peaks in the chromPeaks matrix for each feature).
featureDefinitions(xdata)

# The featureValues method returns a matrix with rows being features and columns
# samples. The content of this matrix can be defined using the value argument.
# The default isvalue = "index" which simply returns the index of the peak (in
# the chromPeaks matrix) assigned to a feature in a sample. Setting value =
# "into" returns a matrix with the integrated signal of the peaks corresponding
# to a feature in a sample. Any column name of the chromPeaks matrix can be
# passed to the argument value. Below we extract the integrated peak intensity
# per feature/sample.

# Extract the into column for each feature.
head(featureValues(xdata, value = "into"))
##       ko15.CDF ko16.CDF ko18.CDF ko19.CDF ko21.CDF ko22.CDF wt15.CDF wt16.CDF
## FT001  1924712  1757151  1714582  1220358  1383417  1180288  2129885  1634342
## FT002   213659   289501   194604    92590       NA   178286   253826   241844
## FT003   349011   451864   337473       NA   343898   208003   364610   360909
## FT004   286221       NA   364300       NA   137261   149098   255698   311297
## FT005   162905       NA   210000       NA   164009   111158       NA       NA
## FT006  1160580       NA  1345515   608017   380970   300716  1286883  1739517
##       wt18.CDF wt19.CDF wt21.CDF wt22.CDF
## FT001  1810112  1507943  1623589  1354005
## FT002   228501   216393   240606   185399
## FT003    54567       NA       NA   221938
## FT004   272268       NA   181640   271128
## FT005       NA       NA   366441       NA
## FT006   515677   621438   639755   508546

# This feature matrix contains NA for samples in which no chromatographic peak
# was detected in the feature’s m/z-rt region. While in many cases there might
# indeed be no peak signal in the respective region, it might also be that there
# is signal, but the peak detection algorithm failed to detect a chromatographic
# peak.

# xcms provides the fillChromPeaks method to fill in intensity data for
# such missing values from the original files. The filled in peaks are added to
# the chromPeaks matrix and are flagged with an 1 in the "is_filled" column.
# Below we perform such a filling-in of missing peaks.
# Filling missing peaks using default settings. Alternatively we could pass a
# FillChromPeaksParam object to the method.
xdata <- fillChromPeaks(xdata)

head(featureValues(xdata))
##       ko15.CDF ko16.CDF ko18.CDF ko19.CDF ko21.CDF ko22.CDF wt15.CDF wt16.CDF
## FT001       68      629     1362     1799     2129     2409     2762     3311
## FT002       48      614     1342     1790     5610     2397     2743     3299
## FT003       30      599     1327     5499     2110     2382     2724     3282
## FT004      107     5368     1398     5500     2140     2425     2816     3329
## FT005      101     5369     1394     5501     2143     2429     5846     5939
## FT006      413     5370     1644     2001     2295     2622     3107     3587
##       wt18.CDF wt19.CDF wt21.CDF wt22.CDF
## FT001     3839     4257     4602     4922
## FT002     3820     4242     4588     4909
## FT003     3796     6112     6227     4891
## FT004     3884     6113     4633     4960
## FT005     6026     6114     4639     6350

# For features without detected peaks in a sample, the method extracts all
# intensities in the mz-rt region of the feature, integrates the signal and adds
# a filled-in peak to the chromPeaks matrix. No peak is added if no signal is
# measured/available for the mz-rt region of the feature. For these, even after
# filling in missing peak data, a NA is reported in the featureValues matrix.

# It should be mentioned that fillChromPeaks uses columns "rtmin", "rtmax",
# "mzmin" and "mzmax" from the featureDefinitions table to define the region
# from which the signal should be integrated to generate the filled-in peak
# signal. These values correspond however to the positions of the peak apex not
# the peak boundaries of all chromatographic peaks assigned to a feature. It
# might be advisable to increase this area in retention time dimension by a
# constant value appropriate to the average peak width in the experiment. Such a
# value can be specified with fixedRt of the FillChromPeaksParam. If the average
# peak width in the experiment is 20 seconds, specifyingfixedRt = 10 ensures
# that the area from which peaks are integrated is at least 20 seconds wide.

# Below we compare the number of missing values before and after filling in
# missing values. We can use the parameter filled of the featureValues method to
# define whether or not filled-in peak values should be returned too.

# Missing values before filling in peaks
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
    FUN = function(z) sum(is.na(z)))
## ko15.CDF ko16.CDF ko18.CDF ko19.CDF ko21.CDF ko22.CDF wt15.CDF wt16.CDF
##       72       74       69      119      142      114      101       96
## wt18.CDF wt19.CDF wt21.CDF wt22.CDF
##       95      128      136      118

# Missing values after filling in peaks
apply(featureValues(xdata), MARGIN = 2, FUN = function(z) sum(is.na(z)))
## ko15.CDF ko16.CDF ko18.CDF ko19.CDF ko21.CDF ko22.CDF wt15.CDF wt16.CDF
##        7        8        4        8       12        8        8        9
## wt18.CDF wt19.CDF wt21.CDF wt22.CDF
##        9       13       13        7

# Next we use the featureSummary function to get a general per-feature summary
# that includes the number of samples in which a peak was found or the number of
# samples in which more than one peak was assigned to the feature. Specifying
# also sample groups breaks down these summary statistics for each individual
# sample group.
head(featureSummary(xdata, group = xdata$sample_group))
##       count      perc multi_count multi_perc       rsd KO_count   KO_perc
## FT001    12 100.00000           0   0.000000 0.1789898        6 100.00000
## FT002    11  91.66667           1   9.090909 0.2409441        5  83.33333
## FT003     9  75.00000           3  33.333333 0.2844117        5  83.33333
## FT004     9  75.00000           2  22.222222 0.3082528        4  66.66667
## FT005     5  41.66667           0   0.000000 0.4824173        4  66.66667
## FT006    11  91.66667           7  63.636364 0.5229604        5  83.33333
##       KO_multi_count KO_multi_perc    KO_rsd WT_count   WT_perc
## FT001              0             0 0.2027357        6 100.00000
## FT002              1            20 0.3653441        6 100.00000
## FT003              2            40 0.2562703        4  66.66667
## FT004              0             0 0.4694612        5  83.33333
## FT005              0             0 0.2492851        1  16.66667
## FT006              4            80 0.5059513        6 100.00000
##       WT_multi_count WT_multi_perc    WT_rsd
## FT001              0             0 0.1601279
## FT002              0             0 0.1069526
## FT003              1            25 0.3335925
## FT004              2            40 0.1840914
## FT005              0             0        NA
## FT006              3            50 0.5758382

# The performance of peak detection, alignment and correspondence should always
# be evaluated by inspecting extracted ion chromatograms e.g. of known
# compounds, internal standards or identified features in general. The
# featureChromatograms function allows to extract chromatograms for each feature
# present in featureDefinitions. The returned Chromatograms object contains an
# ion chromatogram for each feature (each row containing the data for one
# feature) and sample (each column representing containing data for one sample).
# Below we extract the chromatograms for the first 4 features.
feature_chroms <- featureChromatograms(xdata, features = 1:4)
feature_chroms
## Chromatograms with 4 rows and 12 columns
##                   1              2              3              4
##      <Chromatogram> <Chromatogram> <Chromatogram> <Chromatogram>
## [1,]     length: 38     length: 38     length: 39     length: 38
## [2,]     length: 42     length: 41     length: 43     length: 41
## [3,]     length: 36     length: 36     length: 35     length: 35
## [4,]     length: 67     length: 69     length: 70     length: 70
##                   5              6              7              8
##      <Chromatogram> <Chromatogram> <Chromatogram> <Chromatogram>
## [1,]     length: 37     length: 37     length: 39     length: 38
## [2,]     length: 41     length: 41     length: 42     length: 42
## [3,]     length: 36     length: 36     length: 36     length: 36
## [4,]     length: 73     length: 73     length: 66     length: 69
##                   9             10             11             12
##      <Chromatogram> <Chromatogram> <Chromatogram> <Chromatogram>
## [1,]     length: 38     length: 40     length: 38     length: 36
## [2,]     length: 41     length: 44     length: 42     length: 40
## [3,]     length: 37     length: 34     length: 34     length: 36
## [4,]     length: 71     length: 71     length: 75     length: 73
## phenoData with 2 variables
## featureData with 5 variables

# Plot extracted ion chromatograms for feature number 3 and 4
pdf("output/test/eic_features_3and4.pdf")
par(mfrow = c(1, 2))
plot(feature_chroms[3, ], col = group_colors[feature_chroms$sample_group])
plot(feature_chroms[4, ], col = group_colors[feature_chroms$sample_group])
dev.off()

# Finally perform a principal component analysis to evaluate the grouping of the
# samples in this experiment. Note that we did not perform any data
# normalization hence the grouping might (and will) also be influenced by
# technical biases.

# Extract the features and log2 transform them
ft_ints <- log2(featureValues(xdata, value = "into"))

# Perform a PCA omitting all features with an NA in any of the samples. Also,
# the intensities are mean centered.
pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)

## Plot the PCA
cols <- group_colors[xdata$sample_group]
pcSummary <- summary(pc)
pdf("output/test/pca.pdf")
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
    xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
    digits = 3), " % variance"),
    ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
    digits = 3), " % variance"),
    col = "darkgrey", bg = cols, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2], labels = xdata$sample_name, col = "darkgrey",
    pos = 3, cex = 2)
dev.off()

# From the PCA results, we can see the expected separation between the KO and WT
# samples on PC2. On PC1 samples separate based on their ID, samples with an
# ID <= 18 from samples with an ID > 18. This separation might be caused by a
# technical bias (e.g. measurements performed on different days/weeks) or due to
# biological properties of the mice analyzed (sex, age, litter mates etc).
