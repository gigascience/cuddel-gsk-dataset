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
    # peak_features <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_posdata.txt", sep='\t', header=FALSE)
    # # Peak metadata - peak numbers, m/z values, retention times and MPAs
    # peak_meta <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep='\t', header=TRUE)
    # out <- correlatePeaks(peak_features, peak_meta)
    # print(out[1,])
    # write.table(out, file=tmp)
    # tmp_correlatePeaks_output <- read.table("./Study_pos_CorrListP_Routput.txt", sep='\t', header=FALSE)
})

teardown({
    unlink(tmp)
})

test_that("running correlatePeaks produces expected results", {
    # m <- read.table(tmp, sep='\t', header=FALSE)
    # print(tmp_correlatePeaks_output[1,])
    peak_features <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_posdata.txt", sep='\t', header=FALSE)
    # Peak metadata - peak numbers, m/z values, retention times and MPAs
    peak_meta <- read.table("/Users/peterli/PUTMEDID_LCMS_v1.01/Study_pospeaks.txt", sep='\t', header=TRUE)
    out <- correlatePeaks(peak_features, peak_meta)
    # print(toString(out[1,]))

    # Create tab file containg peak comparisons
    write.table(out,
                file = "Study_pos_CorrListP_Routput.txt",
                append = TRUE,
                sep = "\t",
                row.names=FALSE,
                col.names=FALSE,
                quote=FALSE)

    expect_match(toString(out[1,]), "0.74135")
    expect_match(toString(out[19,]), "0.92248")
    expect_match(toString(out[60,]), "0.74666")
    expect_match(toString(out[215,]), "0.72901")
})