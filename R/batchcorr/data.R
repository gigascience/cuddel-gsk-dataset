#' PTnofill
#'
#' Unfilled peaktable from XCMS object containing batch-specific QCs and
#' long-term reference samples from three batches of LC-MS metabolomics data.
#' Peaktable contains missing values, used for the batch alignment algorithm.
#'
#' @format A data frame with 90 observations of 11667 variables
#'   (i.e. metabolite features) - containing missing values
"PTnofill"

#' PTfill
#'
#' Filled peaktable from XCMS object containing batch-specific QCs and
#' long-term reference samples from three batches of LC-MS metabolomics data.
#' Peaktable contains all values, i.e. after hard integration of EICs.
#'
#' @format A data frame with 90 observations of 11667 variables
#'   (i.e. metabolite features) - without missing values
"PTfill"

#' meta
#'
#' Metadata for the injections in PTfill and PTnofill
#'
#' Col 1: Batch identifier
#'
#' Col 2: Sample type identifier (Q for QC and R for Ref)
#'
#' @format A data frame with 90 observations of 2 variables
#'   (i.e. metabolite features) - containing missing values:
#'   \code{batch}, \code{grp}.
"meta"

#' BData
#'
#' Batch object for within batch correction of 'Batch_B'
#'
#' @format A list with QC, Ref and Batch objects
"BData"

#' FData
#'
#' Batch object for within batch correction of 'Batch_F'
#'
#' @format A list with QC, Ref and Batch objects
"FData"

#' HData
#'
#' Batch object for within batch correction of 'Batch_H'
#'
#' @format A list with QC, Ref and Batch objects
"HData"
