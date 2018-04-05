#' Transpose an Assay (Data Frame)
#'
#' @description Transpose an object of class \code{data.frame} that contains
#'    assay measurements while preserving row (feature) and column (sample)
#'    names.
#'
#' @param assay_df A data frame with numeric values to transpose
#' @param firstColIsFeatureNames Are the data row names in the first column of
#'    \code{df}? Defaults to \code{TRUE}. If not, this function assumes that the
#'    data row names are accesible by the \code{\link{rownames}} function.
#'
#' @details This function is designed to transpose "tall" assay data frames
#'    (where genes or proteins are the rows and patient or tumour samples are
#'    the columns). This function also transposes the row (feature) names to
#'    column names and the column (sample) names to row names. Notice that all
#'    rows and columns (other than the feature name column, as applicable) are
#'    numeric.
#'
#' @return The transposition of \code{df}, with row and column names preserved.
#'
#' @export
#'
#' @examples
#'    x_mat <- matrix(rnorm(5000), ncol = 20, nrow = 250)
#'    rownames(x_mat) <- paste0("gene_", 1:250)
#'    colnames(x_mat) <- paste0("sample_", 1:20)
#'    x_df <- as.data.frame(x_mat, row.names = rownames(x_mat))
#'
#'    transpose_assay(x_df, firstColIsFeatureNames = FALSE)
#'
transpose_assay <- function(assay_df, firstColIsFeatureNames = TRUE){

  if(firstColIsFeatureNames){

    featureNames_vec <- assay_df[, 1, drop = TRUE]
    sampleNames_vec <- colnames(assay_df)[-1]

    transpose_df <- as.data.frame(t(assay_df[, -1]))
    rownames(transpose_df) <- NULL

    colnames(transpose_df) <- featureNames_vec

    sampleNames_df <- data.frame(Sample = sampleNames_vec,
                                 stringsAsFactors = FALSE)
    transpose_df <- cbind(sampleNames_df, transpose_df)

  } else {

    featureNames_vec <- rownames(assay_df)
    sampleNames_vec <- colnames(assay_df)

    transpose_df <- as.data.frame(t(assay_df))

    colnames(transpose_df) <- featureNames_vec
    rownames(transpose_df) <- sampleNames_vec

  }

  class(transpose_df) <- class(assay_df)
  transpose_df

}
