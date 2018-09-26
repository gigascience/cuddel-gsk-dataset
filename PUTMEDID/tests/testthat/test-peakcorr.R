# Title     : TODO
# Objective : TODO
# Created by: peterli
# Created on: 10/9/2018

context("running correlatePeaks")

# library(testthat)
# library(PUTMEDID)
# test_check("PUTMEDID")

tmp <- tempfile()
setup({
})

teardown({
    unlink(tmp)
})

test_that("running correlatePeaks produces expected results", {
    peak_features <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_posdata.txt", sep='\t', header=FALSE)
    # Peak metadata - peak numbers, m/z values, retention times and MPAs
    peak_meta <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep='\t', header=TRUE)
    results <- correlatePeaks(peak_features, peak_meta)

    # Create tab file containing peak comparisons
    write.table(results,
                file = tmp,
                append = TRUE,
                sep = "\t",
                row.names=FALSE,
                col.names=FALSE,
                quote=FALSE)

    expect_match(toString(results[1,]), "0.74135")
    expect_match(toString(results[19,]), "0.92248")
    expect_match(toString(results[60,]), "0.74666")
    expect_match(toString(results[215,]), "0.72901")
    expect_match(toString(results[23955,]), "0.94475570")
    expect_match(toString(results[37424,]), "0.9850918")
})