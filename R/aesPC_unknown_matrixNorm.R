#' Normalize and reconstruct the eigenvalues of a data matrix for supervised PCA
#'
#' @description Normalize the columns of a project matrix. For each eigenvector,
#'    swap the signs of the vector elements if the first entry is negative. See
#'    "Details" for more information.
#'
#' @param B A projection matrix: often the matrix of the left singular vectors
#'    given by the Singular Value Decomposition of a data matrix or Grammian.
#' @param d The number of columns of \code{B} to normalize.
#'
#' @return A matrix of the eigenvectors or left singular vectors in \code{B}
#'    transformed to be the left singular values of the original data matrix.
#'
#' @details This function is designed to reconstruct the original first \code{d}
#'    left singular vectors of a data matrix from the first \code{d}
#'    eigenvectors of the Grammian of that data matrix. Basically, after the
#'    data matrix has been centred, the left singular vectors of that data
#'    matrix and the left singular vectors of the Grammian of that data matrix
#'    are equal up to a sign. This function reverses that sign so that the two
#'    sets of singular vectors are equal.
#'
#'    Consider the internal workings of the \code{\link{aespca}} function. This
#'    "sign flipping" changes the eigenvectors of \code{xtx} into the
#'    left singular vectors of \code{scale(X, , center = TRUE, scale = TRUE)}.
#'    Instead of calculating the Grammian, regularising it (by adding some small
#'    \eqn{\lambda} value to the diagonal), taking the SVD of the regularized
#'    Grammian, and extracting the first \eqn{d} eigenvectors, why don't we just
#'    extract the first \eqn{d} singular vectors directly from the scaled data
#'    matrix itself? The regularisation effect only inflates the singular- or
#'    eigen-values anyway, so it has no effect on the singular vectors in any
#'    way. Moreover, the \code{\link{aespca}} function does not even call for
#'    the eigen-values at all, so this whole process is supurfluous. The only
#'    wrinkle is adapting the \code{\link{lars.lsa}} and \code{\link{aespca}}
#'    functions to only operate on the data matrix.
#'
#'    Furthermore, the \code{\link[lars]{lars}} function \emph{can} take in the
#'    full data, instead of just a Grammian. As an enhancement, we should either
#'    update our copy of the lars function in \code{\link{lars.lsa}}, or make a
#'    call to the exported \code{\link[lars]{lars}} function. This will be an
#'    enhancement for a the next version.
#'
#' @seealso \code{\link{aespca}}; \code{\link{lars.lsa}}; \code{\link{AESPCA_pVals}}
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
normalize <- function(B, d){
  # browser()

  # Square root of the column sums of B ^ 2 (the entries squared).
  normB <- sqrt(apply(B ^ 2, 2, sum))
  # Replace any 0 values with 1. Once again, why? The domain is [0, infty), so
  #   shouldn't we replace 0s with the minimum of the non-zero values? Why 1?
  normB[normB == 0] <- 1
  # Divide each column of B by its "norm" value
  B <- t(t(B) / normB)
  # Replace any really small "normed" B values with 0
  B[abs(B) < 1e-6] <- 0

  # This loop moves over columns: it swaps the signs of each value in a column
  #   if the first value in that column is < 0.
  for (i in 1:d){
    # browser()

    # Check for non-"zero" values in the ith column
    tmp <- abs(B[, i]) > 1e-6
    # If we have any non-"zero" entries in the ith column, then
    if(sum(tmp) != 0){

      # we find the indices of these non-small values and
      tmp1 <- B[tmp, i]
      # reverse their sign if the sign of the first non-small entry in B[,i] is
      #   negative.
      B[, i] <- sign(tmp1[1]) * B[, i]

    }

    # END for
  }

  B
}
