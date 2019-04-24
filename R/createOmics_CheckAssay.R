#' Check an Input Assay
#'
#' @description Check the classes, dimensions, missingness, feature variance,
#'   feature type, and feature names of a data frame.
#'
#' @param df An assay data frame supplied to the \code{\link{CreateOmics}}
#'    function. The first column is assumed to be the sample IDs, and will be
#'    ignored. See \code{\link{CheckSampleIDs}} for checking sample IDs.
#' @param removeNear0 Should columns of \code{df} with variance near 0 be
#'    removed? Defaults to \code{TRUE}.
#' @param epsilon Threshold to consider the variance of a column equal to 0.
#'    Defaults to 0.000001.
#'
#' @return The same data frame, without features with 0 variance, \strong{if}
#'    that data frame passes all checks.
#'
#' @details This function checks that the data frame is not a matrix, that the
#'    data frame has more columns than rows (tidy genomic data), that the data
#'    frame contains no missing or character values, that no features of the
#'    data frame have variance less than \code{epsilon} (and removes such
#'    features if \code{removeNear0 = TRUE}), and checks the data frame for
#'    valid column names.
#'
#' @keywords internal
#'
#'
#' @importFrom stats sd
#'
#' @examples
#'  # DO NOT CALL THIS FUNCTION DIRECTLY. CALL FROM WITHIN CreateOmics().
#'
#' \dontrun{
#'  data("colonSurv_df")
#'  CheckAssay(colonSurv_df[, -(1:3)])
#' }
#'
CheckAssay <- function(df, removeNear0 = TRUE, epsilon = 10^-6){
  # browser()


  ###  Check Classes  ###
  if("matrix" %in% class(df) &
     !("data.frame" %in% class(df))){
    stop(
"\n  You have supplied a matrix object to the assayData_df argument. Note that
the pathwayPCA package functions require -Omics data as an N x p data frame
object: this data frame will have one observation per row and one measurement
per column. If your matrix is in 'tall' (p x N) format, please transpose your
matrix with the 't()' function (but pay attention to your column names after
transposition). Next, you can use the 'as.data.frame()' function to transform
your -Omics data matrix to class 'data.frame'. Please see the help information
found in ?CreateOmicsPath for more details. If you have a p x N data frame,
please see the ?TransposeAssay function.")
  }


  ###  Warn for "tall" Data  ###
  warnLvl <- options("warn")$warn
  options(warn = 1)
  if(nrow(df) > ncol(df)){
    warning(
"  The assayData_df argument has more rows than columns. The pathwayPCA
package functions require -Omics data as an N x p data frame object: this data
frame will have one observation per row and one feature per column. If your
assay is in 'tall' (p x N) format, please transpose your assay with the
'TransposeAssay()' function. Please see the ?TransposeAssay function help file
for more information.
")
  }
  options(warn = warnLvl)
  # Warnings are cached and only printed when control returns to the top level.
  # See https://adv-r.hadley.nz/conditions.html for more information.


  ###  Set Aside Sample IDs Column  ###
  outClass <- class(df)
  df2 <- df[, -1]


  ###  Check for Missing or Non-Numeric Values or 0 Variance Features  ###
  if(anyNA(df2)){
    stop("Missing observations are not permitted in the assay data.")
  }
  if(any(vapply(df2, function(x){ !is.numeric(x) }, logical(1)))){
    stop("Non-numeric values are not permitted in the assay data.")
  }

  smallVars <- vapply(df2, sd, numeric(1)) < sqrt(epsilon)
  if(any(smallVars)){

    var0Genes <- colnames(df2)[smallVars]
    message(
      sprintf("%i genes have variance < epsilon and will be removed. These gene(s) are:",
              length(var0Genes))
    )
    print(var0Genes)

    if(removeNear0){
      df2 <- df2[, !smallVars]
    }

  }


  ###  Check Column Names  ###
  bad_names <- .detect_invalid_names(colnames(df2))
  if(length(bad_names) > 0){
    message(
      sprintf("%i gene name(s) are invalid. Invalid name(s) are:",
              length(bad_names))
    )
    print(bad_names)
    message("These genes may be excluded from analysis. Proper gene names
contain alphanumeric characters only, and start with a letter.")
  }


  ###  Return  ###
  out <- cbind(
    df[, 1, drop = FALSE],
    df2
  )
  class(out) <- outClass

  out

}


###  Utility Function  ###
.detect_invalid_names <- function(string_ls){
  # browser()

  unmatchedSingleQ <- grepl("'", string_ls)
  unmatchedDoubleQ <- grepl('"', string_ls)
  nonAlphaNum_idx <- grepl("[^-a-zA-Z0-9-]", string_ls)
  leadingNum_idx <- grepl("^\\d", string_ls)

  bad_idx <- unmatchedSingleQ +
    unmatchedDoubleQ +
    nonAlphaNum_idx +
    leadingNum_idx
  bad_idx <- as.logical(bad_idx)

  string_ls[bad_idx]

}
