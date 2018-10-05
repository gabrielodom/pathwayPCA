#' Load Test Data onto Training-Data Estimated PCs
#'
#' @description Given a list of loading vectors from a training data set,
#'    calculate the PCs of the test data set.
#'
#' @param design_df A test data frame with rows as samples and named features
#'    as columns
#' @param loadings_ls A list of \eqn{p \times d} loading vectors or matrices as
#'    returned by either the \code{\link{SuperPCA_pVals}},
#'    \code{\link{AESPCA_pVals}}, or \code{\link{ExtractAESPCs}} functions.
#'    These lists of loadings will have feature names as their row names. Such
#'    feature names must match a subset of the column names of \code{design_df}
#'    exactly, as pathway-specific test-data subsetting is performed by column
#'    name.
#'
#' @return A data frame with the PCs from each pathway concatenated by column.
#'    If you have the \code{tidyverse} loaded, this object will display as a
#'    \code{\link[tibble]{tibble}}.
#'
#' @details This function takes in a list of loadings and a training-centered
#'    test data set, applies over the list of loadings, subsets the columns of
#'    the test data by the row names of the loading vectors, right-multiplies
#'    the test-data subset matrix by the loading vector / matrix, and returns a
#'    data frame of the test-data PCs for each loading vector.
#'
#' @export
#'
#' @examples
#'    NULL
LoadOntoPCs <- function(design_df, loadings_ls){
  # browser()

  ###  PCs List  ###
  PCs_ls <- lapply(loadings_ls, function(x){

    tryCatch(
      {
        genes <- rownames(x)
        testSubset <- as.matrix(design_df[genes])
        testSubset %*% x
      },
      error = function(e) NULL
    )

  })

  ###  Check for Missing PCs  ###
  error_idx <- sapply(PCs_ls, is.null)
  if(any(error_idx)){

    message(sprintf("PC extraction failed for %i pathway(s). These pathway(s) are:",
                    sum(error_idx)))
    print(which(error_idx))

    out_mat <- do.call(cbind, PCs_ls[!error_idx])
    names_char <- names(PCs_ls)[!error_idx]

  } else {

    out_mat <- do.call(cbind, PCs_ls)
    names_char <- names(PCs_ls)

  }

  ###  Return  ###
  out_df <- as.data.frame(out_mat)
  colnames(out_df) <- names_char
  class(out_df) <- c("tbl_df", "tbl", "data.frame")
  out_df
}
