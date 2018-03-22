#' Normalize a matrix for supervised PCA
#'
#' @param B A matrix.
#' @param d The number of columns of \code{B} to normalize.
#'
#' @description A function that norms a matrix, but I don't understand any of
#'   it.
#'
#' @return A modified version of the \code{B} matrix.
#'
#' @details I met with James and Steven on 26 September and neither of them
#'   understood the sign reversal in the last line of the internal \code{for()}
#'   loop. Based on how it's called in the \code{\link{aespca}} function, it has
#'   something to do with adjusting the AES-PCA eigenvectors returned by the
#'   \code{\link{lars.lsa}} function. DOCUMENT THIS.
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
normalize <- function(B, d){
  # browser()

  # Square root of the column sums of B ^ 2 (the entries squared). Why though?
  normB <- sqrt(apply(B ^ 2, 2, sum))
  # Replace any 0 values with 1. Once again, why? The domain is [0, infty), so
  #   shouldn't we replace 0s with the minimum of the non-zero values? Why 1?
  normB[normB == 0] <- 1
  # Divide each column of B by its "norm" value
  B <- t(t(B) / normB)
  # Replace any really small "normed" B values with 0
  B[abs(B) < 1e-6] <- 0

  # This loop moves over columns: it swaps the signs of each value in a column
  #   if the first value in that column is < 0. I still don't get what the point
  #   is of what we're doing.
  for (i in 1:d){
    # browser()

    # Check for non-"zero" values in the ith column
    tmp <- abs(B[, i]) > 1e-6
    # If we have any non-"zero" entries in the ith column, then
    if(sum(tmp) != 0){

      # we find the indices of these non-small values and
      tmp1 <- B[tmp, i]
      # reverse their sign if the sign of the first non-small entry in B[,i] is
      #   negative. What is going on?? This is to negate the sign swap from the
      #   Grammian matrix. Now that we have changed to the svd of the data, we
      #   should be able to remove it.
      B[, i] <- sign(tmp1[1]) * B[, i]

    }

    # END for
  }

  B
}
