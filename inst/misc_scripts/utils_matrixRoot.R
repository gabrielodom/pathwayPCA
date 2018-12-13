#' Positive root of a symmetric matrix
#'
#' @description Calculate the matrix root of a symmetric matrix via the Spectral
#'   Decomposition.
#'
#' @param x A symmetric (necessarily square) matrix.
#' @param root A positive real number.
#'
#' @return A matrix, that when multiplied by itself \code{root} times, yields
#'   \code{x}.
#'
#' @details This function decomposes \code{x} into \eqn{V \times D \times V^T}
#'   via the \code{\link[base]{eigen}} function, sets any numerically negative
#'   eigenvalues to 0, calculates the root of these eigenvalues as \eqn{D^r},
#'   then returns the matrix \eqn{V \times D^r \times V^T}.
#'
#'   See \url{https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix}.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   X <- matrix(rnorm(25), ncol = 5);   xTx <- t(X) %*% X
#'   matrixRoot(xTx)
#'   matrixRoot(xTx, root = 3)
matrixRoot <- function(x, root = 2){

  # x.eigen <- eigen(x)
  # d <- x.eigen$values
  # # Zero-out any negative eigenvalues
  # d <- (d + abs(d))/2
  # v <- x.eigen$vectors
  #
  # if(root == 2){
  #   return(v %*% diag(sqrt(d)) %*% t(v))
  # } else {
  #   return(v %*% diag(d ^ (1 / root)) %*% t(v))
  # }

}
