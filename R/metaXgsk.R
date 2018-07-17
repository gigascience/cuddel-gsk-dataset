library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(mice)
library(ggfortify)
library(metaX)
library(VIM)
library(ropls)

# Set working directory
wd <- setwd("")

# Let R know where GSK data set is located
datadir <- "/home/peter/"

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
neg_peaks <- neg_peaklist[ , !names(neg_peaklist) %in% c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "block1neg", "block2neg", "block3neg", "block4neg", "percent_nas", "isotopes", "adduct", "pcgroup")]
# Rename idx column to name
colnames(neg_peaks)[1] <- "name"

#############################
# Load peak data into metaX #
#############################

# Output peak intensity data from xcms object
write.table(neg_peaks, file = "new_peaks.csv", sep =",", row.names = TRUE, col.names = TRUE)

# Create new metaXpara object
para <- new("metaXpara")
# Set output directory
outdir(para) <- "test"
prefix(para) <- "gsk_"

# Import data from XCMS into metaX
para <- importDataFromXCMS(para, file="neg_peaks.csv")
# Check data import
head(para@rawPeaks[,1:20])
name GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_12r GSK_neg_block1_13r GSK_neg_block1_14r GSK_neg_block1_15r GSK_neg_block1_16r
1    1           314657.3           253010.0           234015.0           229210.3           271279.5           225345.9           290000.1
2    2          2214161.1          2213167.3          2262108.8          2285202.3          1992895.3          2245666.6          2200922.3

###########################
# Create sample list file #
###########################

# Sample file needs to look like this:
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
class <- as.character(class)
order <- meta_all[, "order"]
sample_meta <- cbind(sample, batch, class, order)
sample_meta <- as.matrix(sample_meta)

write.table(sample_meta, file = "sample_meta.tab", sep ="\t", row.names = FALSE, col.names = TRUE)
sfile <- paste(wd, "/sample_meta.tab", sep="")
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)


######################################################
# Remove metabolite features detected in <50% of QCs #
######################################################

para <- filterQCPeaks(para, ratio = 0.5)


#########################
# Impute missing values #
#########################

para <- missingValueImpute(para, method = "knn")


###################
# Remove outliers #
###################

# para <- autoRemoveOutlier(para, valueID="value")


##################
# Normalise data #
##################

para <- doQCRLSC(para, cpu = 3)

# Check normalisation
getPeaksTable(para, valueID="value")
getPeaksTable(para, valueID="valueNorm")

#######
# PCA #
#######

para <- transformation(para$metaXpara, valueID = "valueNorm")
metaX::plotPCA(para, valueID = "valueNorm", scale = "pareto", center = TRUE, rmQC = FALSE)

##########
# PLS-DA #
##########

para <- transformation(para$metaXpara, valueID = "valueNorm")
plsdaPara <- new("plsDAPara")
plsdaPara@nperm <- 100
plsda.res <- runPLSDA(para = para, plsdaPara = plsdaPara, sample = c("NA", "A", "B", "C", "D"), valueID = "valueNorm")

#########
# ROPLS #
#########

norm_peaks <- getPeaksTable(para, valueID="valueNorm")
# Create dataMatrix object
dataMatrix <- norm_peaks[c(1,5:ncol(norm_peaks))]
rownames(dataMatrix) <- dataMatrix[,1]
dataMatrix <- dataMatrix[ -c(1) ]

# Create sampleMetadata object which needs to look like this:
## > head(sampleMetadata)
## age   bmi gender
## HU_011  29 19.75      M
## HU_014  59 22.64      F
## HU_015  42 22.72      M
## HU_017  41 23.03      M
## HU_018  34 20.96      M
## HU_019  35 23.41      M

sampleMetadata <- norm_peaks[c(2:4)]
rownames(sampleMetadata) <- norm_peaks[,1]

# variableMetadata object contains the chemical identity of the peaks
# Have not created this variableMetadata object

# Run PCA
gsk.pca <- opls(dataMatrix)

# Save PCA plot
pdf('gsk_pca.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("overview", "outlier", "x-score", "x-loading"))
plot(gsk.pca, typeVc = typeC, parDevNewL = FALSE)
dev.off()

# Check if there is an partitioning based on regimen class
regimenFc <- sampleMetadata[, "class"]
pdf('gsk_pca_regimen.pdf')
plot(gsk.pca, typeVc = "x-score", parAsColFcVn = regimenFc, parEllipsesL = TRUE, parDevNewL = FALSE)
dev.off()

# Perform PLS and PLS-DA
gsk.plsda <- opls(dataMatrix, regimenFc, permI = 100, plotL = FALSE)
pdf('gsk_plsda_regimen.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("permutation", "overview", "outlier", "x-score"))
plot(gsk.plsda, typeVc = typeC, parDevNewL = FALSE)
dev.off()

#########################################
# metaXpipe whole data analysis process #
#########################################

p <- metaXpipe(para, plsdaPara=plsdaPara, cvFilter=0.3, remveOutlier = TRUE, outTol = 1.2, doQA = TRUE, doROC = TRUE, qcsc = 1, nor.method = "pqn", pclean = TRUE, t = 1, scale = "uv",)


