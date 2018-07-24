# Title     : PCA and CVA analysis on GSK data
# Objective : Perform PCA and CVA analysis on GSK data using Morpho R package
# Created by: peterli
# Created on: 23/7/2018

library(xcms)
library(CAMERA)
library(naturalsort)
library(matrixStats)
library(mice)
library(ggfortify)
library(metaX)
library(VIM)
library(ropls)
library(MAIT)
library(Morpho)

# Working directory is /home/peter/metax/gsk/neg/morpho
wd <- setwd("/home/peter/metax/gsk/neg/morpho")
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

file_name_neg <- meta_all[,"file_name_neg"]
file_name_neg <- as.character(file_name_neg)

# Create new metaXpara object
para <- new("metaXpara")
# Set output directory
outdir(para) <- "test"
prefix(para) <- "gsk"

para <- importDataFromXCMS(para, file="xcms_neg_peaks.csv")
# Check data import
head(para@rawPeaks[,1:20])
## name GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_12r GSK_neg_block1_13r GSK_neg_block1_14r GSK_neg_block1_15r GSK_neg_block1_16r
## 1    1           314657.3           253010.0           234015.0           229210.3           271279.5           225345.9           290000.1
## 2    2          2214161.1          2213167.3          2262108.8          2285202.3          1992895.3          2245666.6          2200922.3

###########################
# Create sample list file #
###########################

# Sample file needs to look like this:
## sample  batch   class   order
## batch01_QC01    1       NA      1
## batch01_QC02    1       NA      2
## batch01_QC03    1       NA      3
## batch01_C05     1       C       4
## batch01_S07     1       S       5
## batch01_C10     1       C       6

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

#################################################
# Do PCA on unprocessed negative QC and samples #
#################################################

xcms_neg_peaks = read.csv("xcms_neg_peaks.csv")

# Prepare neg_peaklist data
pca_data <- xcms_neg_peaks[, file_name_neg]

# Prepare block information for labelling data points
block <- integer(0)
sample_names <- colnames(pca_data)
for (i in 1:length(sample_names)) {
    if (grepl("block1", sample_names[i]) == 1) {
        block[i] <- "block1"
    }
    else if (grepl("block2", sample_names[i]) == 1) {
        block[i] <- "block2"
    }
    else if (grepl("block3", sample_names[i]) == 1) {
        block[i] <- "block3"
    }
    else if (grepl("block4", sample_names[i]) == 1) {
        block[i] <- "block4"
    }
}

# Prepare QC and sample information for labelling data points
pca_meta_qc_sample <- integer(0)
for (i in 1:length(sample_names)) {
    if (meta_all[which(meta_all[,"file_name_neg"]==sample_names[i]), "type"] == 'QC') {
        pca_meta_qc_sample[i] <- "QC"
    }
    else {
        pca_meta_qc_sample[i] <- "Sample"
    }
}

# Transpose data
pca_data <- t(pca_data)
# Add block information, QC and sample names to PCA data
pca_data <- cbind(pca_data, block, pca_meta_qc_sample)
# PCA cannot be performed on data with missing values
# Therefore remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:409]), data = pca_data, shape= "block", colour = 'pca_meta_qc_sample', main = 'PCA on unprocessed negative QC and sample data')
ggsave("unprocessed_neg_qc_sample_pca.png")


######################################################
# Remove metabolite features detected in <50% of QCs #
######################################################

para <- filterQCPeaks(para, ratio = 0.5)
## Remove peaks which the percent is more than 0.5 with intensity are NA!
## 1216
## Save the removed peaks to file: test/gsk-filterQCPeaks


###########################
# Do PCA on filtered data #
###########################

# Data for PCA analysis needs to be like this:
head(pca_data)
##                          X2      X11      X12      X46      X47      X51      X52      X55      X67      X69      X84      X85      X92      X94
## GSK_neg_block1_09r 20.99499 15.88976 19.82944 16.27867 24.90300 15.27127 24.30450 21.24730 18.92648 17.00591 21.09528 19.81756 17.83334 16.91873
## GSK_neg_block1_10r 20.99518 16.06161 19.77672 16.28347 24.88628 15.14840 24.29889 21.17415 19.13544 16.80561 21.10195 19.80083 17.72526 16.79727
## GSK_neg_block1_16r 20.99132 15.93550 19.71773 16.31846 24.91606 15.47256 24.31950 21.27124 19.02973 16.75619 21.02317 19.85418 17.85304 17.05119
## GSK_neg_block1_21r 21.01415 15.77040 19.70531 15.99287 24.79537 14.89362 24.20281 21.07380 19.05643 16.64718 20.95879 19.80568 17.50214 16.86835
## GSK_neg_block1_26r 21.03203 16.05183 19.83050 16.42897 24.95190 15.51925 24.35288 21.30907 19.10666 16.79962 21.17225 19.81461 17.81682 16.87366
## GSK_neg_block1_31r 20.93778 15.59914 19.72883 15.83856 24.73029 14.58911 24.12010 21.19329 18.95975 16.72093 20.93511 19.81138 17.60189 16.73402
##                       X2983    X2986  block pca_meta_qc_sample regimens
## GSK_neg_block1_09r 15.97102 16.53889 block1                 QC       QC
## GSK_neg_block1_10r 15.97276 16.25647 block1                 QC       QC
## GSK_neg_block1_16r 15.76988 16.26237 block1                 QC       QC
## GSK_neg_block1_21r 15.99016 16.50822 block1                 QC       QC
## GSK_neg_block1_26r 15.71526 16.18567 block1                 QC       QC
## GSK_neg_block1_31r 15.96108 16.17907 block1                 QC       QC

# Prepare data
mydata <- getPeaksTable(para, valueID="value")
mydata[1:3, 1:6]
##               sample class batch order        1      10
## 1 GSK_neg_block1_09r    QC     1     1 314657.3 2796343
## 2 GSK_neg_block1_10r    QC     1     2 253010.0 2448407
## 3 GSK_neg_block1_12r     C     1     3 234015.0 2865853

# Use sample names for column names
rownames(mydata) <- mydata[,"sample"]
# Remove not required columns
nr <- c("sample", "class", "batch", "order")
nr <- colnames(mydata) %in% c("sample", "class", "batch", "order")
pca_data <- mydata[, !nr]

# Prepare block information for labelling data points
block <- integer(0)
sample_names <- rownames(pca_data)
for (i in 1:length(sample_names)) {
    if (grepl("block1", sample_names[i]) == 1) {
        block[i] <- "block1"
    }
    else if (grepl("block2", sample_names[i]) == 1) {
        block[i] <- "block2"
    }
    else if (grepl("block3", sample_names[i]) == 1) {
        block[i] <- "block3"
    }
    else if (grepl("block4", sample_names[i]) == 1) {
        block[i] <- "block4"
    }
}

# Transpose data
#pca_data <- t(pca_data)
# Add block information to PCA data
pca_data <- cbind(pca_data, block)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour= "block", main = 'PCA on filtered data before normalisation')
ggsave("filtered_before_norm_pca.png")


#########################
# Impute missing values #
#########################

para <- missingValueImpute(para, method = "knn")
## missingValueImpute: value
## Wed Jul 18 14:12:08 2018        Missing value imputation for 'value'
## Missing value in total: 118587
## Missing value in QC sample: 21461
## Missing value in non-QC sample: 97126
## Wed Jul 18 14:12:08 2018        The ratio of missing value: 16.6654%
## Cluster size 1831 broken into 1821 10
## Cluster size 1821 broken into 1800 21
## Cluster size 1800 broken into 87 1713
## Done cluster 87
## Cluster size 1713 broken into 1467 246
## Done cluster 1467
## Done cluster 246
## Done cluster 1713
## Done cluster 1800
## Done cluster 21
## Done cluster 1821
## Done cluster 10
## <=0: 142
## Missing value in total after missing value inputation: 0
## <=0 value in total after missing value inputation: 0
## Warning message:
## In knnimp(x, k, maxmiss = rowmax, maxp = maxp) :
## 87 rows with more than 50 % entries missing;
## mean imputation used for these rows



###################
# Remove outliers #
###################

# para <- autoRemoveOutlier(para, valueID="value")


##################
# Normalise data #
##################

para <- doQCRLSC(para, cpu = 12)

para <- transformation(para$metaXpara, valueID = "valueNorm")

# Check normalisation
getPeaksTable(para, valueID="value")
getPeaksTable(para, valueID="valueNorm")


#################################
# Do PCA on non-normalised data #
#################################

# Data for PCA analysis need to be like this:
head(pca_data)
##                          X2      X11      X12      X46      X47      X51      X52      X55      X67      X69      X84      X85      X92      X94
## GSK_neg_block1_09r 20.99499 15.88976 19.82944 16.27867 24.90300 15.27127 24.30450 21.24730 18.92648 17.00591 21.09528 19.81756 17.83334 16.91873
## GSK_neg_block1_10r 20.99518 16.06161 19.77672 16.28347 24.88628 15.14840 24.29889 21.17415 19.13544 16.80561 21.10195 19.80083 17.72526 16.79727
## GSK_neg_block1_16r 20.99132 15.93550 19.71773 16.31846 24.91606 15.47256 24.31950 21.27124 19.02973 16.75619 21.02317 19.85418 17.85304 17.05119
## GSK_neg_block1_21r 21.01415 15.77040 19.70531 15.99287 24.79537 14.89362 24.20281 21.07380 19.05643 16.64718 20.95879 19.80568 17.50214 16.86835
## GSK_neg_block1_26r 21.03203 16.05183 19.83050 16.42897 24.95190 15.51925 24.35288 21.30907 19.10666 16.79962 21.17225 19.81461 17.81682 16.87366
## GSK_neg_block1_31r 20.93778 15.59914 19.72883 15.83856 24.73029 14.58911 24.12010 21.19329 18.95975 16.72093 20.93511 19.81138 17.60189 16.73402
##                       X2983    X2986  block pca_meta_qc_sample regimens
## GSK_neg_block1_09r 15.97102 16.53889 block1                 QC       QC
## GSK_neg_block1_10r 15.97276 16.25647 block1                 QC       QC
## GSK_neg_block1_16r 15.76988 16.26237 block1                 QC       QC
## GSK_neg_block1_21r 15.99016 16.50822 block1                 QC       QC
## GSK_neg_block1_26r 15.71526 16.18567 block1                 QC       QC
## GSK_neg_block1_31r 15.96108 16.17907 block1                 QC       QC

# Prepare data
mydata <- getPeaksTable(para, valueID="value")
mydata[1:3,1:6]
##               sample class batch order        1      10
## 1 GSK_neg_block1_09r    QC     1     1 314657.3 2796343
## 2 GSK_neg_block1_10r    QC     1     2 253010.0 2448407
## 3 GSK_neg_block1_12r     C     1     3 234015.0 2865853

# Use sample names for column names
rownames(mydata) <- mydata[,"sample"]
# Remove not required columns
nr <- c("sample", "class", "batch", "order")
nr <- colnames(mydata) %in% c("sample", "class", "batch", "order")
pca_data <- mydata[, !nr]

# Prepare block information for labelling data points
block <- integer(0)
sample_names <- rownames(pca_data)
for (i in 1:length(sample_names)) {
    if (grepl("block1", sample_names[i]) == 1) {
        block[i] <- "block1"
    }
    else if (grepl("block2", sample_names[i]) == 1) {
        block[i] <- "block2"
    }
    else if (grepl("block3", sample_names[i]) == 1) {
        block[i] <- "block3"
    }
    else if (grepl("block4", sample_names[i]) == 1) {
        block[i] <- "block4"
    }
}

# Transpose data
#pca_data <- t(pca_data)
# Add block information to PCA data
pca_data <- cbind(pca_data, block)
# PCA cannot be performed on data with missing values
# Remove peak rows if it contains missing values
#pca_data <- pca_data[ , colSums(is.na(pca_data)) == 0]

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour= "block", main = 'PCA before normalisation')
ggsave("before_norm_pca.png")

#####################################
# Do PCA to check signal correction #
#####################################

> head(rsd_knn_qc_neg[,1:4])
GSK_neg_block1_09r GSK_neg_block1_10r GSK_neg_block1_16r GSK_neg_block1_21r
2          2214161.07         2213167.30         2200922.28         2229604.08
11           80404.29           90185.92           80860.06           70636.49
12         1052155.89         1013478.90          968455.80          955786.38
46           95204.15           95296.60           96498.84           76137.69
47        37436017.50        36930624.08        37326859.84        34001468.40
51           47954.21           53991.18           53033.27           34404.22

> head(peak_data)
ID             sample       value batch class order
1  2 GSK_neg_block1_09r  2214161.07     1    NA     1
2 11 GSK_neg_block1_09r    80404.29     1    NA     1
3 12 GSK_neg_block1_09r  1052155.89     1    NA     1
4 46 GSK_neg_block1_09r    95204.15     1    NA     1
5 47 GSK_neg_block1_09r 37436017.50     1    NA     1
6 51 GSK_neg_block1_09r    47954.21     1    NA     1

# Get signal corrected peak data
para <- transformation(para$metaXpara, valueID = "valueNorm")
norm_data <- getPeaksTable(para, valueID = "valueNorm")
# Use sample column as rownames
rownames(norm_data) <- norm_data$sample
# Remove sample, class, batch and order columns
drop_cols <- c("sample", "class", "batch", "order")
norm_data <- norm_data[ , !(names(norm_data) %in% drop_cols)]
# Transpose data
norm_data <- t(norm_data)
# Order data based on rownames (which are peak IDs)
norm_data <- norm_data[naturalorder(rownames(norm_data)), ]

pca_data <- norm_data
# Prepare block information for labelling data points
block <- integer(0)
sample_names <- colnames(pca_data)
for (i in 1:length(sample_names)) {
  if (grepl("block1", sample_names[i]) == 1) {
    block[i] <- "block1"
  }
  else if (grepl("block2", sample_names[i]) == 1) {
    block[i] <- "block2"
  }
  else if (grepl("block3", sample_names[i]) == 1) {
    block[i] <- "block3"
  }
  else if (grepl("block4", sample_names[i]) == 1) {
    block[i] <- "block4"
  }
}

# Prepare QC and sample information for labelling data points
pca_meta_qc_sample <- integer(0)
for (i in 1:length(sample_names)) {
  if (meta_all[which(meta_all[,"file_name_neg"]==sample_names[i]), "type"] == 'QC') {
    pca_meta_qc_sample[i] <- "QC"
  }
  else {
    pca_meta_qc_sample[i] <- "Sample"
  }
}

regimens <- lapply(sample_names, function(sample_name) {
  regimen <- meta_all[which(meta_all[,"file_name_neg"] == sample_name), "Regimen"]
  if(is.na(regimen)) {
    regimen <- "QC"
  }
    regimen <- as.character(regimen)
  })
regimens <- unlist(regimens)

# Transpose data
pca_data <- t(pca_data)
# Add block information, QC and sample names to PCA data
pca_data <- cbind(pca_data, block, pca_meta_qc_sample, regimens)

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")

pca_results <- prcomp(pca_data[,1:506])

# Highlight sample type and batch in PCA plot
autoplot(prcomp(pca_data[,1:506]), data = pca_data, shape= "block", colour = 'pca_meta_qc_sample', main = 'PCA on normalised negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("norm_neg_qc_sample_pca.png")

# Highlight regimens in PCA plot
autoplot(prcomp(pca_data[,1:506]), data = pca_data, colour = 'regimens', main = 'PCA on normalised negative QC and sample data', frame = TRUE, frame.type = 'norm')
ggsave("norm_neg_qc_sample_pca.png")




#######
# CVA #
#######

cvall <- CVA(pca_results$x[, 1:5], groups)
head(proc$PCscores[, 1:5])
## PC1          PC2          PC3          PC4          PC5
## skull_0144_ch_fe 0.032385579 -0.006687871 -0.019071891 -0.010079025 -0.020413868
## skull_0149_ch_fe 0.053375114 -0.010885558 -0.013503190  0.009357859 -0.005599557
## skull_0150_ch_fe 0.047636357 -0.028018847 -0.004376455  0.007873201  0.016918763
## skull_0152_ch_fe 0.081384812 -0.034833964  0.031245020  0.006269748  0.015588755
## skull_0157_ch_fe 0.004453751 -0.014276104 -0.001606326  0.016094292 -0.012025064
## skull_0164_ch_fe 0.074665565  0.008102766  0.009828817  0.009840437 -0.001317569


#######
# PCA #
#######

data(iris)
vari <- iris[, 1:4]
head(vari)
## Sepal.Length Sepal.Width Petal.Length Petal.Width
## 1          5.1         3.5          1.4         0.2
## 2          4.9         3.0          1.4         0.2
## 3          4.7         3.2          1.3         0.2
## 4          4.6         3.1          1.5         0.2
## 5          5.0         3.6          1.4         0.2
nrow(vari)
## [1] 150

facto <- iris[, 5]
head(facto)
## [1] setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica
length(facto)
## [1] 150
pca.1 <-groupPCA(vari, groups=facto, rounds=0, mc.cores=12)


para <- transformation(para$metaXpara, valueID = "valueNorm")
#vari <- getPeaksTable(para, valueID="valueNorm")
vari <- getPeaksTable(para, valueID="value")
#vari <- vari[1:150, 5:ncol(vari)]
vari <- vari[1:150, 5:8]
vari <- signif(vari, digits=2)
vari <- as.matrix(vari)
vari[1:5,1:5]
##          1      10     1003     1004     1005
## 1 314657.3 2796343 287626.5 491842.1 68558.91
## 2 253010.0 2448407 290437.1 483219.8 62892.78
## 3 234015.0 2865853 321979.1 394396.6 57336.35
## 4 229210.3 2737243 318245.3 351762.9 81769.64
## 5 271279.5 2819283 201149.5 390935.0 58849.65

nrow(vari)
## [1] 371

facto <- sample_meta[, "batch"]
facto <- as.factor(facto)
length(facto)
## [1] 371

head(facto)
## [1] 1 1 1 1 1 1
## Levels: 1 2 3 4


pca.1 <- groupPCA(vari, groups=facto, rounds=0, mc.cores=12)

pdf('scatterplotMatrix.pdf')
scatterplotMatrix(pca.1$Scores, groups=facto, ellipse=TRUE, by.groups=TRUE, var.labels=c("PC1", "PC2", "PC3"))
dev.off()


#######
# CVA #
#######

groups <- as.factor(regimens)
cvall <- CVA(pca_results$x[, 1:5], groups)

# plot the CVA
pdf('regimens-cva.pdf')
plot(cvall$CVscores, col=facto, pch=as.numeric(facto), typ="n", asp=1,
xlab=paste("1st canonical axis", paste(round(cvall$Var[1, 2], 1), "%")),
ylab=paste("2nd canonical axis", paste(round(cvall$Var[2, 2], 1), "%")))
text(cvall$CVscores, as.character(facto), col=as.numeric(facto), cex=.7)
dev.off()

cvall <- CVA(proc$PCscores[, 1:5], groups)
head(proc$PCscores[, 1:5])
## PC1          PC2          PC3          PC4          PC5
## skull_0144_ch_fe 0.032385579 -0.006687871 -0.019071891 -0.010079025 -0.020413868
## skull_0149_ch_fe 0.053375114 -0.010885558 -0.013503190  0.009357859 -0.005599557
## skull_0150_ch_fe 0.047636357 -0.028018847 -0.004376455  0.007873201  0.016918763
## skull_0152_ch_fe 0.081384812 -0.034833964  0.031245020  0.006269748  0.015588755
## skull_0157_ch_fe 0.004453751 -0.014276104 -0.001606326  0.016094292 -0.012025064
## skull_0164_ch_fe 0.074665565  0.008102766  0.009828817  0.009840437 -0.001317569


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
dataMatrix <- norm_peaks[c(1, 5:ncol(norm_peaks))]
rownames(dataMatrix) <- dataMatrix[, 1]
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
rownames(sampleMetadata) <- norm_peaks[, 1]

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


###########################################
# Perform univariate statistical analysis #
###########################################

# Confirm regimens in excel spreadsheet
## A = food yes, exercise no
## B = food no, exercise no
## C = food yes, exercise yes
## D = food no, exercise yes

# However, in PhD thesis it is:
## A = food yes, exercise no
## B = food no, exercise yes
## C = food yes, exercise yes
## D = food no, exercise no

ratioPairs(para) <- "A:B"
# Check norm values are in para object
getPeaksTable(para, valueID="valueNorm")

#addValueNorm(para) <- para
para <- transformation(para$metaXpara, valueID = "valueNorm")
plsdaPara <- new("plsDAPara")
plsdaPara@nperm <- 10
plsda.res <- runPLSDA(para = para, plsdaPara = plsdaPara, sample = c("NA", "A", "B", "C", "D"), valueID = "valueNorm")
res <- peakStat(para, plsdaPara, doROC = TRUE)

#########################################
# metaXpipe whole data analysis process #
#########################################

p <- metaXpipe(para, plsdaPara=plsdaPara, cvFilter=0.3, remveOutlier = TRUE, outTol = 1.2, doQA = TRUE, doROC = TRUE, qcsc = 1, nor.method = "pqn", pclean = TRUE, t = 1, scale = "uv",)



