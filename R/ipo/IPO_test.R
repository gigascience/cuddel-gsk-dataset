# Title     : IPO_test.R
# Objective : Test use of IPO
# Created by: peterli
# Created on: 1/11/2018

library(IPO)

datapath <- system.file("cdf", package = "faahKO")
datafiles <- list.files(datapath, recursive = TRUE, full.names = TRUE)

# Need to choose between centWave or matchedFilter
peakpickingParameters <- getDefaultXcmsSetStartingParams('matchedFilter')
# Setting levels for step to 0.2 and 0.3 (hence 0.25 is the center point)
peakpickingParameters$step <- c(0.2, 0.3)
peakpickingParameters$fwhm <- c(40, 50)
# Setting only one value for steps therefore this parameter is not optimized
peakpickingParameters$steps <- 2

# Optimise peak picking parameters
time.xcmsSet <- system.time({ # measuring time
    resultPeakpicking <-
        optimizeXcmsSet(files = datafiles[1:2],
        params = peakpickingParameters,
        BPPARAM=bpparam(),
        nSlaves = 0,  # nSlaves set to 0 to stop IPO nSlaves and xcms BPPARAM conflict
        subdir = NULL,
        plot = TRUE)
})

# Optimize retention time correction and grouping parameters
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
retcorGroupParameters$gapExtend <- 2.7
time.RetGroup <- system.time({ # measuring time
    resultRetcorGroup <- optimizeRetGroup(xset = optimizedXcmsSetObject,
        params = retcorGroupParameters,
        BPPARAM=bpparam(),
        nSlaves = 0,
        subdir = NULL,
        plot = TRUE)
    })
