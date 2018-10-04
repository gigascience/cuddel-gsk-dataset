# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 26/9/2018

context("running annotateMass")

# library(testthat)
# library(PUTMEDID)
# test_check("PUTMEDID")

tmp <- tempfile()
setup({
})

teardown({
    unlink(tmp)
})

test_that("running annotateMass has correct parameters", {
    # peak_features <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_posdata.txt", sep='\t', header=FALSE)
    # # Peak metadata - peak numbers, m/z values, retention times and MPAs
    # peak_meta <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep='\t', header=TRUE)
    results <- annotateMassmatch()

    # Create tab file containing peak comparisons
    # write.table(results,
    #     file = "tmp.txt",
    #     append = TRUE,
    #     sep = "\t",
    #     row.names=FALSE,
    #     col.names=FALSE,
    #     quote=FALSE)

    expect_match(toString(results[1,]), "pos")      # elabel
    expect_match(toString(results[2,]), "3")        # lowval
    expect_match(toString(results[3,]), "0.00075")  # limtol
    expect_match(toString(results[4,]), "2.5")      # rtlim
    expect_match(toString(results[5,]), "30")       # rtmin
    expect_match(toString(results[6,]), "1200")     # rtman
    expect_match(toString(results[7,]), "0.9")      # mycorrlim
    expect_match(toString(results[8,]), "1")        # etol
    expect_match(toString(results[9,]), "E")        # PeakLabel[3269]
    expect_match(toString(results[10,]), "A")       # PeakLabel[3270]
    expect_match(toString(results[11,]), "3802")    # ordpname[3802]
    expect_match(toString(results[12,]), "3802")    # sortPeaks[3802]
    expect_match(toString(results[13,]), "0.15815") # final tdiff
    expect_match(toString(results[14,]), "30.01056")  # MFdata[1]
    ncorr <- results[15,]
    expect_match(ncorr, "23888")  # ncorr
    totCorrs_size <- results[16,]
    expect_match(totCorrs_size, "37424")  # totCorrs_size
    totCorrs_1 <- results[17,]
    expect_match(totCorrs_1, "50457")  # totCorrs_1
    finCorrs_size <- results[18,]
    expect_match(finCorrs_size, "23888")  # finCorrs_size
    finCorrs_2 <- results[19,]
    expect_match(finCorrs_2, "50465")  # finCorrs_2
    ogno <- results[20,]
    expect_match(ogno, "306")  # ogno
    massmatch_3802 <- results[21,]
    expect_match(massmatch_3802, "993.6790501")  # massmatch_3802
    Inum_3801 <- results[22,]
    expect_match(Inum_3801, "19.4998498891375")  # Inum_3801
    Amode_15 <- results[23,]
    expect_match(Amode_15, "5")  # Amode_15
})
