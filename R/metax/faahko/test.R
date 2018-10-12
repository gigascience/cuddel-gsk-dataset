library("metaX")

list.files(system.file("cdf", package = "faahKO"), recursive = TRUE, full.names = TRUE)
[1] "/usr/lib64/R/library/faahKO/cdf/KO/ko15.CDF"
[2] "/usr/lib64/R/library/faahKO/cdf/KO/ko16.CDF"
[3] "/usr/lib64/R/library/faahKO/cdf/KO/ko18.CDF"
[4] "/usr/lib64/R/library/faahKO/cdf/KO/ko19.CDF"
[5] "/usr/lib64/R/library/faahKO/cdf/KO/ko21.CDF"
[6] "/usr/lib64/R/library/faahKO/cdf/KO/ko22.CDF"
[7] "/usr/lib64/R/library/faahKO/cdf/WT/wt15.CDF"
[8] "/usr/lib64/R/library/faahKO/cdf/WT/wt16.CDF"
[9] "/usr/lib64/R/library/faahKO/cdf/WT/wt18.CDF"
[10] "/usr/lib64/R/library/faahKO/cdf/WT/wt19.CDF"
[11] "/usr/lib64/R/library/faahKO/cdf/WT/wt21.CDF"
[12] "/usr/lib64/R/library/faahKO/cdf/WT/wt22.CDF"

# Use metaXpara class object to manage file path info and parameters for data processing
para <- new("metaXpara")
## set the MS data path
dir.case(para) <- system.file("cdf/KO", package = "faahKO")
dir.ctrl(para) <- system.file("cdf/WT", package = "faahKO")

# Set sample list file path
sampleFile <- system.file("extdata/faahKO_sampleList.txt", package = "metaX")
sampleListFile(para) <- sampleFile
samList <- read.delim(sampleFile)
print(samList)

# Set params for peak peaking
xcmsSet.peakwidth(para) <- c(20,50)
xcmsSet.snthresh(para) <- 10
xcmsSet.prefilter(para) <- c(3,100)
xcmsSet.noise(para) <- 0
xcmsSet.nSlaves(para) <- 4

# Set params for peak annotation
group.bw(para) <- 5
group.minfrac(para) <- 0.3
group.mzwid(para) <- 0.015
group.max(para) <- 1000

# Do peak picking and annotation using peakFinder() function
p <- peakFinder(para)