packages1 <- c("Nozzle.R1", "ggplot2", "reshape2", "plyr", "BBmisc", "mixOmics", "missForest", "doParallel", "DiscriMiner", "ape")
packages2 <- c("scatterplot3d", "pheatmap", "bootstrap", "boot", "caret", "dplyr", "stringr", "RColorBrewer", "RCurl", "lattice")
packages3 <- c("lattice", "data.table", "tidyr", "scales", "VennDiagram", "pROC", "readr", "e1071", "randomForest", "coop")
packages4 <- c("fpc", "ROCR", "devtools", "DiffCorr", "XML", "RSQLite", "igraph")

bioc_packages1 <- c("pcaMethods", "multtest", "vsn", "pls", "mzID", "MSnbase", "faahKO", "mzR")
bioc_packages2 <- c("xcms", "preprocessCore", "impute", "CAMERA", "ropls", "sva", "SSPA", "rgl")

remove.packages(packages1)
remove.packages(packages2)
remove.packages(packages3)
remove.packages(packages4)

remove.packages(bioc_packages1)
remove.packages(bioc_packages2)

remove.packages("metaX")
