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

# Create new metaXpara object
para <- new("metaXpara")
# Set output directory
outdir(para) <- "test"
prefix(para) <- "gsk"

# Import data from XCMS into metaX
para <- importDataFromXCMS(para, file="neg_peaks.csv")
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


######################################################
# Remove metabolite features detected in <50% of QCs #
######################################################

para <- filterQCPeaks(para, ratio = 0.5)
## Remove peaks which the percent is more than 0.5 with intensity are NA!
## 1216
## Save the removed peaks to file: test/gsk-filterQCPeaks


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

para <- doQCRLSC(para, cpu = 8)

# Check normalisation
getPeaksTable(para, valueID="value")
getPeaksTable(para, valueID="valueNorm")


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



