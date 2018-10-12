library(ropls)

# Load sacurine data
data(sacurine)

# Attach dataset to search path
attach(sacurine)

# Display summary of data using strF function from ropls package
strF(dataMatrix)
## dim  class    mode typeof   size NAs  min mean median max
## 183 x 109 matrix numeric double 0.2 Mb   0 -0.3  4.2    4.3   6
## (2-methoxyethoxy)propanoic acid isomer (gamma)Glu-Leu/Ile ... Valerylglycine isomer 2  Xanthosine
## HU_011                            3.019766011        3.888479324 ...             3.889078716 4.075879575
## HU_014                             3.81433889        4.277148905 ...             4.181765852 4.195761901
## ...                                       ...                ... ...                     ...         ...
## HU_208                            3.748127215        4.523763202 ...             4.634338821 4.487781609
## HU_209                            4.208859398        4.675880567 ...              4.47194762 4.222953354

strF(sampleMetadata)
## age     bmi gender
## numeric numeric factor
## nRow nCol size NAs
## 183    3 0 Mb   0
## age   bmi gender
## HU_011   29 19.75      M
## HU_014   59 22.64      F
## ...     ...   ...    ...
## HU_208   27 18.61      F
## HU_209 17.5 21.48      F

strF(variableMetadata)
## msiLevel      hmdb chemicalClass
## numeric character     character
## nRow nCol size NAs
## 109    3 0 Mb   0
## msiLevel      hmdb chemicalClass
## (2-methoxyethoxy)propanoic acid isomer        2                  Organi
## (gamma)Glu-Leu/Ile                            2                  AA-pep
## ...                                         ...       ...           ...
## Valerylglycine isomer 2                       2           AA-pep:AcyGly
## Xanthosine                                    1 HMDB00299        Nucleo

# Perform PCA
sacurine.pca <- opls(dataMatrix)
## PCA
## 183 samples x 109 variables
## standard scaling of predictors
##       R2X(cum) pre ort
## Total    0.501   8   0

# Save PCA summary plot as PDF
pdf('pca.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("overview", "outlier", "x-score", "x-loading"))
plot(sacurine.pca, typeVc = typeC, parDevNewL = FALSE)
dev.off()

# Check if there is an partitioning based on gender
genderFc <- sampleMetadata[, "gender"]
pdf('pca_gender.pdf')
plot(sacurine.pca, typeVc = "x-score", parAsColFcVn = genderFc, parEllipsesL = TRUE, parDevNewL = FALSE)
dev.off()

# Perform PLS and PLS-DA
sacurine.plsda <- opls(dataMatrix, genderFc, plotL = FALSE)
pdf('plsda_gender.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("permutation", "overview", "outlier", "x-score"))
plot(sacurine.plsda, typeVc = typeC, parDevNewL = FALSE)
dev.off()

# Perform OPLS and OPLS-DA
sacurine.oplsda <- opls(dataMatrix, genderFc, predI = 1, orthoI = NA, plotL = FALSE)
pdf('oplsda_gender.pdf')
layout(matrix(1:4, nrow = 2, byrow = TRUE))
for(typeC in c("permutation", "overview", "outlier", "x-score"))
plot(sacurine.oplsda, typeVc = typeC, parDevNewL = FALSE)
dev.off()