library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(mice)
library(ggfortify)
library(metaX)
library(VIM)

# Let R know where GSK data set is located
datadir = "/home/peter/"

######################################
# Read in metadata for data analysis #
######################################

# Read metadata files
meta <- read.csv(paste(datadir, "gsk/meta/meta.csv", sep =""))
meta_all <- read.csv(paste(datadir, "gsk/meta/meta_all.csv", sep =""))
# Sort by file name neg
meta <- meta[naturalorder(meta$file_name_neg),]
meta_all <- meta_all[naturalorder(meta_all$file_name_neg),]

# Create file paths for negative files
neg_dir = paste(datadir, "gsk/raw/esi_neg/netcdf", sep ="")
file_name_neg <- meta_all[,"file_name_neg"]
file_name_neg <- as.character(file_name_neg)
neg_file_paths <- lapply(file_name_neg, function (x) {
    if (grepl("block1", x) == 1) {
        return(paste(neg_dir, "/block1neg/", x, ".cdf", sep = ""))
    }
    else if (grepl("block2", x) == 1) {
        return(paste(neg_dir, "/block2neg/", x, ".cdf", sep = ""))
    }
    else if (grepl("block3", x) == 1) {
        return(paste(neg_dir, "/block3neg/", x, ".cdf", sep = ""))
    }
    else if (grepl("block4", x) == 1) {
        return(paste(neg_dir, "/block4neg/", x, ".cdf", sep = ""))
    }
})
neg_file_paths <- as.character(neg_file_paths)


######################################
# Apply XCMS onto QC and sample data #
######################################

# Create xcmsSet object using findPeaks parameters from Eva's thesis
# This is compute intensive and takes time to complete!!
neg_xset <- xcmsSet(neg_file_paths, step = 0.02, snthresh=3, mzdiff = 0.05)

# Match peaks representing same analyte across samples
grp_neg_xset <- group(neg_xset, bw = 10, mzwid = 0.05)

# Create CAMERA object
neg_xsa <- xsAnnotate(grp_neg_xset)
# Get peaklist
neg_peaklist <- getPeaklist(neg_xsa)
# Add numeric index to negative peaklist to keep track of peaks during data pre-treatment
neg_peaklist$idx <- seq.int(nrow(neg_peaklist))
# Move index to left hand side of data frame
neg_peaklist <- neg_peaklist[, c(ncol(neg_peaklist), 1:(ncol(neg_peaklist)-1))]

# The peak list needs to look like this from the metaX tutorial:
head(para@rawPeaks[,1:4])
name batch01_QC01 batch01_QC02 batch01_QC03
1  78.02055 14023.0      13071.0     15270.0
2 452.00345 22455.0      10737.0     27397.0
3 138.96337 6635.4       8062.3      6294.6

# However, the peak list currently looks like this:
> head(neg_peaklist[,1:20])
idx       mz    mzmin    mzmax        rt     rtmin     rtmax npeaks block1neg block2neg block3neg block4neg percent_nas GSK_neg_block1_09r
1   1 57.97605 57.97542 57.97685   64.9456   59.1596   92.6520    369        92        86        93        92    2.083333           314657.3
2   2 59.01426 59.01381 59.01453   64.6192   62.7960   67.1400    371        92        93        93        93    0.000000          2214161.1
3   3 61.98876 61.98832 61.98905 1393.2574 1390.5089 1395.7195    371        92        93        93        93    0.000000          4573570.5

# List of columns to delete
# del_cols <- c("idx", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks")
rawPeaks <- neg_peaklist[ , !names(neg_peaklist) %in% c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "block1neg", "block2neg", "block3neg", "block4neg", "percent_nas")]
# Rename idx column to name
colnames(rawPeaks)[1] <- "name"

#############################
# Load peak data into metaX #
#############################

# Output peak intensity data from xcms object
write.table(rawPeaks, file = "rawPeaks.tab", sep ="\t", row.names = TRUE, col.names = TRUE)

# Create new metaXpara object
para <- new("metaXpara")

# Import data from XCMS into metaX
para <- importDataFromXCMS(para, file="rawPeaks.csv")
# Check data import
head(para@rawPeaks[,1:20])
idx       mz    mzmin    mzmax        rt     rtmin     rtmax npeaks block1neg block2neg block3neg block4neg percent_nas GSK_neg_block1_09r
1   1 57.97605 57.97542 57.97685   64.9456   59.1596   92.6520    369        92        86        93        92    2.083333           314657.3
2   2 59.01426 59.01381 59.01453   64.6192   62.7960   67.1400    371        92        93        93        93    0.000000          2214161.1

# Create sample list file that looks like this:
sample  batch   class   order
batch01_QC01    1       NA      1
batch01_QC02    1       NA      2
batch01_QC03    1       NA      3
batch01_C05     1       C       4
batch01_S07     1       S       5
batch01_C10     1       C       6

# Get batch metadata from sample names
sample <- meta_all[, "file_name_neg"]
sample <- as.character(sample)
batch <- meta_all[, "block"]
class <- meta_all[, "Regimen"]
order <- meta_all[, "order"]
sampleListFile <- cbind(sample, batch, class, order)
sampleListFile <- as.matrix(sampleListFile)

write.table(sampleListFile, file = "sampleListFile.tab", sep ="\t", row.names = FALSE, col.names = TRUE)
sfile <- paste(datadir, "gsk/raw/esi_neg/netcdf/sampleListFile.tab", sep="")
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)

######################################################
# Remove metabolite features detected in <50% of QCs #
######################################################

para <- filterQCPeaks(para, ratio = 0.5)

#########################
# Impute missing values #
#########################

para <- missingValueImpute(para,method="knn")


##################
# Normalise data #
##################

res <- doQCRLSC(para, cpu=1)



