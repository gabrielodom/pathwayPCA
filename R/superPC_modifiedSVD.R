#' Function Title
#'
#' @param mat A matrix of data frame in "tall" format (p x n)
#' @param n.components How many singular values / vectors to return? Must be an
#'   integer less than min(p, n). Best performance increase is for values much
#'   less than min(p, n). Defaults to \code{NULL}.
#'
#' @description Center and compute the fast SVD of the matrix \code{mat}
#'
#' @return A list containing:
#' \itemize{
#'   \item{u : }{The first \code{n.components} left singular vectors of
#'     \code{mat}.}
#'   \item{d : }{The largest \code{n.component} singular vectors of \code{mat}.}
#'   \item{v : }{The first \code{n.components} right singular vectors of
#'     \code{mat}.}
#'   \item{feature.means : }{A named vector of the feature means of \code{mat}.}
#' }
#'
#' @details The \code{mysvd} function takes in a tall -omics data matrix,
#'   extracts the feature means, centers the matrix on this mean vector, and
#'   calculates the Singular Value Decomposition of the centered data matrix.
#'   Currently, the SVD is calculated via the \code{\link[corpcor]{fast.svd}}
#'   function from \code{corpcor} package. However, this function calculates
#'   all the singular vectors, even when \code{n.components} is non-\code{NULL}.
#'   We should update this to the \code{\link[rsvd]{rsvd}} function from the
#'   \code{rsvd} package.
#'
#' @export
#'
#' @importFrom corpcor fast.svd
#'
#' @examples
#'   NULL

mysvd <- function(mat, n.components = NULL){

  # finds PCs of matrix x
  p <- nrow(mat)
  n <- ncol(mat)

  # center the observations (rows)
  feature.means <- rowMeans(mat)
  mat <- t(scale(t(mat),
                 center = feature.means,
                 scale = FALSE))


  if(is.null(n.components)){
    n.components <- min(n, p)
  }

  junk <- corpcor::fast.svd(mat)
  nc <- min(ncol(junk$u), n.components)

  return(list(u = junk$u[, 1:nc],
              d = junk$d[1:nc],
              v = junk$v[, 1:nc],
              feature.means = feature.means))

}
