# Title     : List_CorrData.R
# Objective : Calculation of correlation between peaks
# Created by: peterli
# Created on: 10/9/2018

library('pracma')

correlatePeaks <- function(peakFeatures, mass_data, calccorr="P", corrcheck="0.7", rtread="5") {
    # Check if peakFeatures has correct data structure corresponding to mass_data
    dimensionality <- nrow(peakFeatures)
    #print(paste0("dimensionality = ", dimensionality))
    ncols <- ncol(peakFeatures)
    #print(paste0("ncols = ", ncols))
    if (calccorr == "P" | calccorr == "p" | calccorr == "S" | calccorr == "s") {
        # Data needs to be n rows (peaks) x ncols columns (samples)
        # In beanshell peakFeatures has header so has extra row, it is ((ncols+1 == dimensionality) || (dimensionality+1 != npeaks))
        if (ncols+1 == dimensionality-1 | dimensionality != nrow(mass_data)) {
            stop("Incorrect data to calculate correlation - program aborted")
        }
    } else {
        # In beanshell peakFeatures has header so has extra row, it is (ncols != dimensionality) || (dimensionality+1 != npeaks)
        if (ncols != dimensionality-1 | dimensionality != nrow(mass_data)) {
            stop("Data should be a square matrix (n x n) with n = no. of peaks - program aborted")
        }
    }

    # Set up internal variables
    rtdiff <- as.double(rtread)
    corrlim <- as.double(corrcheck)
    tcorr <- 0.0
    tempStr <- ""
    OKdata <- TRUE
    npeaks <- nrow(mass_data)
    # print(paste0("npeaks = ", npeaks))

    # To hold retention times
    rt <- vector(mode="double", length=npeaks)
    rt <- mass_data[, 3]
    # To hold peak numbers
    PeakNo <- vector(mode="integer", length=npeaks)
    PeakNo <- mass_data[, 1]
    # To hold peak numbers as double objects
    PeakNoDble <- vector(mode="double", length=npeaks)
    PeakNoDble <- as.double(PeakNo)

    # Sort peak numbers in order
    order <- shellSort(PeakNoDble)

    # To hold correlation results between 2 peaks
    listofdata <- mat.or.vec(0, 3)

    # Transpose peak data
    values <- t(peakFeatures)

    rtgap <- 0
    # Select a peak
    for(i in 1:npeaks) {
        j <- i + 1
        # print(paste0("i: ", i))
        # Select another peak to compare with first selected peak
        for(j in j:npeaks) {
            if(j > npeaks)
                break
            # print(paste0("j: ", j))
            tcorr <- 0
            # Uses Peak No. order to calculate difference in RT
            rtgap=rt[order[i]]-rt[order[j]]
            if (rtgap<0) {
                rtgap <- -rtgap
            }
            # Check if RT gap is more than the RT difference threshold
            if (rtgap<rtdiff) {
                # Do Pearson correlation calculation if more that RT threshold
                if (calccorr == "P" | calccorr == "p") {
                    corr <- cor.test(values[, order[i]], values[, order[j]], method = "pearson")
                    tcorr <- corr$estimate
                }
                # Do Spearman correlation calculation if more than RT threshold
                else if (calccorr == "S" | calccorr == "s") {
                    corr <- cor.test(values[, order[i]], values[, order[j]], method = "spearman")
                    tcorr <- corr$estimate
                }
                else {
                    tcorr <- values[i][j]
                }
            }

            # If correlation coefficient is more than correlation threshold
            # then include it in results
            if (tcorr>corrlim) {
                if (PeakNo[order[i]]<PeakNo[order[j]]) {
                    row <- c(PeakNo[order[i]], PeakNo[order[j]], tcorr)
                    listofdata <- rbind(listofdata, row)
                }
                else {
                    row <- c(PeakNo[order[j]], PeakNo[order[i]], tcorr)
                    listofdata <- rbind(listofdata, row)
                }
            }
        }
    }

    corrlist <- listofdata
    return(corrlist)
}
