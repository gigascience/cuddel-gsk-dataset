# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 18/7/2018

library(Morpho)
library(shapes)

## all examples are kindly provided by Marta Rufino

if (require(shapes)) {
    # perform procrustes fit on raw data
    alldat <- procSym(abind(gorf.dat, gorm.dat))
    # create factors
    groups <- as.factor(c(rep("female", 30),rep("male", 29)))
    # perform CVA and test Mahalanobis distance
    # between groups with permutation test by 100 rounds)
    cvall <- CVA(alldat$orpdata, groups, rounds=10000)
    ## visualize a shape change from score -5 to 5:
    cvvis5 <- 5*matrix(cvall$CVvis[, 1], nrow(cvall$Grandm), ncol(cvall$Grandm)) + cvall$Grandm
    cvvisNeg5 <- -5*matrix(cvall$CVvis[, 1], nrow(cvall$Grandm), ncol(cvall$Grandm)) + cvall$Grandm
    plot(cvvis5, asp=1)
    points(cvvisNeg5, col=2)
    for (i in 1:nrow(cvvisNeg5))
    lines(rbind(cvvis5[i, ], cvvisNeg5[i, ]))
}
### Morpho CVA
data(iris)
vari <- iris[, 1:4]
facto <- iris[, 5]

cva.1 = CVA(vari, groups=facto)
## get the typicality probabilities and resulting classifications - tagging
## all specimens with a probability of < 0.01 as outliers (assigned to no class)
typprobs <- typprobClass(cva.1$CVscores, groups=facto)
print(typprobs)
## visualize the CV scores by their groups estimated from (cross-validated)
## typicality probabilities:
if (require(car)) {
    scatterplot(cva.1$CVscores[, 1], cva.1$CVscores[, 2], groups=typprobs$groupaffinCV, smooth=FALSE, reg.line=FALSE)
}
# plot the CVA
pdf('iris-cva.pdf')
plot(cva.1$CVscores, col=facto, pch=as.numeric(facto), typ="n", asp=1,
    xlab=paste("1st canonical axis", paste(round(cva.1$Var[1, 2], 1), "%")),
    ylab=paste("2nd canonical axis", paste(round(cva.1$Var[2, 2], 1), "%")))
text(cva.1$CVscores, as.character(facto), col=as.numeric(facto), cex=.7)
dev.off()

# add chull (merge groups)
for(jj in 1:length(levels(facto))){
    ii = levels(facto)[jj]
    kk = chull(cva.1$CVscores[facto == ii, 1:2])
    lines(cva.1$CVscores[facto == ii,1][c(kk, kk[1])], cva.1$CVscores[facto == ii, 2][c(kk, kk[1])], col=jj)
}

# add 80% ellipses
if (require(car)) {
    for(ii in 1:length(levels(facto))){
        dataEllipse(cva.1$CVscores[facto == levels(facto)[ii], 1],
        cva.1$CVscores[facto == levels(facto)[ii], 2],
        add = TRUE, levels = .80, col = c(1:7)[ii])}
}

# histogram per group
if (require(lattice)) {
    histogram(~cva.1$CVscores[, 1]|facto,
    layout = c(1, length(levels(facto))),
    xlab = paste("1st canonical axis", paste(round(cva.1$Var[1, 2], 1),"%")))
    histogram(~cva.1$CVscores[, 2]|facto, layout = c(1,length(levels(facto))),
    xlab = paste("2nd canonical axis", paste(round(cva.1$Var[2, 2], 1), "%")))
}

# plot Mahalahobis
dendroS = hclust(cva.1$Dist$GroupdistMaha)
dendroS$labels = levels(facto)
par(mar = c(4, 4.5, 1, 1))
dendroS = as.dendrogram(dendroS)
pdf('iris-cva-dendro.pdf')
plot(dendroS, main = '', sub = '', xlab = "Geographic areas", ylab = 'Mahalahobis distance')
dev.off()

# Variance explained by the canonical roots:
cva.1$Var
# or plot it:
pdf('iris-cva-barplot.pdf')
barplot(cva.1$Var[, 2])
dev.off()

# another landmark based example in 3D:
data(boneData)
groups <- name2factor(boneLM, which = 3:4)
proc <- procSym(boneLM)
cvall <- CVA(proc$orpdata, groups)
#' ## visualize a shape change from score -5 to 5:
cvvis5 <- 5*matrix(cvall$CVvis[, 1], nrow(cvall$Grandm), ncol(cvall$Grandm)) + cvall$Grandm
cvvisNeg5 <- -5*matrix(cvall$CVvis[, 1], nrow(cvall$Grandm), ncol(cvall$Grandm)) + cvall$Grandm
## Not run:
# visualize it
deformGrid3d(cvvis5, cvvisNeg5, ngrid = 0)

## End(Not run)

#for using (e.g. the first 5) PCscores, one will do:
cvall <- CVA(proc$PCscores[, 1:5], groups)
#' ## visualize a shape change from score -5 to 5:
cvvis5 <- 5*cvall$CVvis[, 1] + cvall$Grandm
cvvisNeg5 <- -5*cvall$CVvis[, 1] + cvall$Grandm
cvvis5 <- showPC(cvvis5, proc$PCs[, 1:5], proc$mshape)
cvvisNeg5 <- showPC(cvvisNeg5, proc$PCs[, 1:5], proc$mshape)
## Not run:
#visualize it
deformGrid3d(cvvis5, cvvisNeg5, ngrid = 0)

## End(Not run)