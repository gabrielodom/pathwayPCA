#' Positive Root of a Symmetric Matrix
#'
#' @description Calculate the matrix root of a symmetric matrix via the Spectral
#'   Decomposition.
#'
#' @param x A symmetric (necessarily square) matrix
#' @param root A positive number
#'
#' @return A matrix, that when multiplied by itself \code{root} times, yields
#'   \code{x}.
#'
#' @details This function decomposes \code{x} into $V x D x V^T$ via the
#'   \code{eigen()} function, sets any numerically negative eigenvalues to 0,
#'   calculates the root of these eigenvalues as D^*, then returns the matrix
#'   $V x D^* x V^T$.
#'
#' @export
#'
#' @examples
#'   X <- matrix(rnorm(25), ncol = 5);   xTx <- t(X) %*% X
#'   matrixRoot(xTx)
#'   matrixRoot(xTx, root = 3)
matrixRoot <- function(x, root = 2){

  x.eigen <- eigen(x)
  d <- x.eigen$values
  # Zero-out any negative eigenvalues
  d <- (d + abs(d))/2
  v <- x.eigen$vectors

  if(root == 2){
    return(v %*% diag(sqrt(d)) %*% t(v))
  } else {
    return(v %*% diag(d ^ (1 / root)) %*% t(v))
  }

}
