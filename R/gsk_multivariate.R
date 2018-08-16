# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 16/8/2018

library(Morpho)

#############################################################################
# Perform Canonical Variate Analysis on food and exercise regimens with QCs #
#############################################################################

ht(iris)
## Sepal.Length Sepal.Width Petal.Length Petal.Width   Species
## 1            5.1         3.5          1.4         0.2    setosa
## 2            4.9         3.0          1.4         0.2    setosa
## 3            4.7         3.2          1.3         0.2    setosa
## 4            4.6         3.1          1.5         0.2    setosa
## 5            5.0         3.6          1.4         0.2    setosa
## 6            5.4         3.9          1.7         0.4    setosa
## 7            4.6         3.4          1.4         0.3    setosa
## 8            5.0         3.4          1.5         0.2    setosa
## 9            4.4         2.9          1.4         0.2    setosa
## 10           4.9         3.1          1.5         0.1    setosa
## 141          6.7         3.1          5.6         2.4 virginica
## 142          6.9         3.1          5.1         2.3 virginica
## 143          5.8         2.7          5.1         1.9 virginica
## 144          6.8         3.2          5.9         2.3 virginica
## 145          6.7         3.3          5.7         2.5 virginica
## 146          6.7         3.0          5.2         2.3 virginica
## 147          6.3         2.5          5.0         1.9 virginica
## 148          6.5         3.0          5.2         2.0 virginica
## 149          6.2         3.4          5.4         2.3 virginica
## 150          5.9         3.0          5.1         1.8 virginica

vari <- iris[,1:4]
h(vari)
## Sepal.Length Sepal.Width Petal.Length Petal.Width
## 1          5.1         3.5          1.4         0.2
## 2          4.9         3.0          1.4         0.2
## 3          4.7         3.2          1.3         0.2
## 4          4.6         3.1          1.5         0.2
## 5          5.0         3.6          1.4         0.2
## 6          5.4         3.9          1.7         0.4

facto <- iris[,5]
h(facto)
## [1] setosa setosa setosa setosa setosa setosa
## Levels: setosa versicolor virginica

pca_results <- prcomp(pca_data[, 1:506])

# Do CVA
groups <- as.factor(regimens)
cvall <- CVA(pca_results$x[, 1:5], groups)

# Plot CVA results
pdf('regimens-cva.pdf')
#png(file = "myplot.png", bg = "white")
plot(cvall$CVscores, col=groups, pch=as.numeric(groups), typ="n", asp=1,
xlab=paste("1st canonical axis", paste(round(cvall$Var[1, 2], 1), "%")),
ylab=paste("2nd canonical axis", paste(round(cvall$Var[2, 2], 1), "%")))
text(cvall$CVscores, as.character(groups), col=as.numeric(groups), cex=.7)

# add chull (merge groups)
for(jj in 1:length(levels(groups))){
    ii=levels(groups)[jj]
    kk=chull(cvall$CVscores[groups==ii,1:2])
    lines(cvall$CVscores[groups==ii,1][c(kk, kk[1])],
    cvall$CVscores[groups==ii, 2][c(kk, kk[1])], col=jj)
}

# add 80% ellipses
if (require(car)) {
    for(ii in 1:length(levels(groups))){
        dataEllipse(cvall$CVscores[groups==levels(groups)[ii],1],
        cvall$CVscores[groups==levels(groups)[ii],2],
        add=TRUE, levels=.80, col=c(1:7)[ii])}
}

dev.off()


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
ropls.norm.data.pca <- opls(dataMatrix)

# Save PCA plot
pdf('gsk_pca.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("overview", "outlier", "x-score", "x-loading"))
plot(ropls.norm.data.pca, typeVc = typeC, parDevNewL = FALSE)
dev.off()

# Check if there is an partitioning based on regimen class
regimenFc <- sampleMetadata[, "class"]
pdf('gsk_pca_regimen.pdf')
plot(ropls.norm.data.pca, typeVc = "x-score", parAsColFcVn = regimenFc, parEllipsesL = TRUE, parDevNewL = FALSE)
dev.off()

# Perform PLS and PLS-DA
gsk.plsda <- opls(dataMatrix, regimenFc, permI = 100, plotL = FALSE)
pdf('gsk_plsda_regimen.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("permutation", "overview", "outlier", "x-score"))
plot(gsk.plsda, typeVc = typeC, parDevNewL = FALSE)
dev.off()


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


#########################
# Normalise data to QCs #
#########################

# Get QC data
qc_names <- meta_all[meta_all[, "type"] == "QC", "file_name_neg"]
qc_names <- as.character(qc_names)
qc_data <- sig_corr_qc_sample_neg[ ,qc_names]

# Calculate average value for each QC peak
qc_means <- rowMeans(qc_data)

# Divide all data with QC data
qcnorm_qc_sample_neg <- sig_corr_qc_sample_neg[ , qc_sample_names]
# new_pca_data <- matrix(nrow = nrow(pca_data), ncol = ncol(pca_data))
# colnames(new_pca_data) <- colnames(pca_data)
# rownames(new_pca_data) <- rownames(pca_data)

# Test
qcnorm_qc_sample_neg <- cbind(qcnorm_qc_sample_neg, qc_means)
qcnorm_qc_sample_neg <- as.matrix(qcnorm_qc_sample_neg)
counter <- 0
ans <- apply(qcnorm_qc_sample_neg, 1, function(x) {
    counter <<- counter + 1
    # print(x)
    # print(paste("Length of x:", length(x)))
    qc_average <- as.numeric(x[length(x)])
    # print("##### End of matrix row #####")
    # for(i in 1:length(x)-1) {
    # print(x[length(x)])
    #   print(x[i] / x[length(x)])
    # }
    stuff <- lapply(x, function (y) {
        # print(paste("y:", y))
        # qc_average <- as.numeric(y[length(y)])
        # print(paste("qc_average: ", qc_average))
        new_value <- y / qc_average
        # print(paste("new value: ", new_value))
        # print("#####")
        return(new_value)
    })
    # print(paste("stuff: ", stuff))
    # print(stuff)
    stuff <- as.numeric(stuff)
    # print(paste("length of stuff: ", length(stuff)))
    # print(class(stuff))
    # print(paste("counter:", counter))
    qcnorm_qc_sample_neg[counter,] <<- stuff
})

# Clean up
qcnorm_qc_sample_neg <- qcnorm_qc_sample_neg[, -ncol(qcnorm_qc_sample_neg)]

# Do PCA
pca_data <- qcnorm_qc_sample_neg

# Transpose data
pca_data <- t(pca_data)
# Add regimen information to PCA data
pca_data <- cbind(pca_data, regimen)

write.table(pca_data, file = "pca_data.csv", sep =",", row.names = TRUE, col.names = TRUE)
pca_data <- read.table(file = "pca_data.csv", sep=",")
autoplot(prcomp(pca_data[,1:ncol(pca_data)-1]), data = pca_data, colour = 'regimen', main = 'PCA on QC-normalised negative QC and sample data')
ggsave("pca_qc_normalised.png")

