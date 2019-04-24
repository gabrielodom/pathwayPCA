#' Singular Value Decomposition wrapper for supervised PCA
#'
#' @param mat A matrix of data frame in "tall" format (\eqn{p \times n}).
#' @param method What function should be used to extract the left- and right-
#'   singular vectors and singular values? Any function that returns the values
#'   as a list with components \code{u}, \code{v}, and \code{d} is appropriate.
#'   Defaults to \code{\link[base]{svd}}.
#' @param n.components How many singular values / vectors to return? Must be an
#'   integer less than \eqn{min(p, n)}. Best performance increase is for values
#'   much less than \eqn{min(p, n)}. Defaults to \code{NULL}.
#'
#' @description Center and compute the SVD of a matrix
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{u} : }{The first \code{n.components} left singular vectors of
#'      \code{mat}.}
#'   \item{\code{d} : }{The largest \code{n.component} singular values of
#'      \code{mat}.}
#'   \item{\code{v} : }{The first \code{n.components} right singular vectors of
#'      \code{mat}.}
#'   \item{\code{feature.means} : }{A named vector of the feature means of
#'      \code{mat}.}
#' }
#'
#' @details The \code{mysvd} function takes in a tall -Omics data matrix,
#'   extracts the feature means, centers the matrix on this mean vector, and
#'   calculates the Singular Value Decomposition (SVD) of the centered data
#'   matrix. Currently, the SVD is calculated via the
#'   \code{\link[corpcor]{fast.svd}} function from \code{corpcor} package.
#'   However, this function calculates all the singular vectors, even when
#'   \code{n.components} is non-\code{NULL}. We should experiment with other SVD
#'   functions, such as the \code{\link[rsvd]{rsvd}} function from the
#'   \code{rsvd} package. ENHANCEMENT.
#'
#' @keywords internal
#'
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use SuperPCA_pVals() instead
#'   
#' \dontrun{
#'   data("colon_pathwayCollection")
#'   data("colonSurv_df")
#'   
#'   colon_OmicsSurv <- CreateOmics(
#'     assayData_df = colonSurv_df[,-(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "surv"
#'   )
#'   
#'   asthmaGenes_char <-
#'     getTrimPathwayCollection(colon_OmicsSurv)[["KEGG_ASTHMA"]]$IDs
#'   
#'   mysvd(t(getAssay(colon_OmicsSurv))[asthmaGenes_char, ])
#' }
#'   

mysvd <- function(mat, method = svd, n.components = NULL){

  # finds PCs of matrix x
  p <- nrow(mat)
  n <- ncol(mat)

  # center the observations (rows)
  feature.means <- rowMeans(mat)
  # mat <- t(scale(t(mat), center = feature.means, scale = FALSE))
  # We moved the scaling step to object creation. See CreateOmics()
  mat <- as.matrix(mat)


  if(is.null(n.components)){
    n.components <- min(n, p)
  }

  junk <- method(mat)
  nc <- min(ncol(junk$u), n.components)

  list(
    u = junk$u[, seq_len(nc)],
    d = junk$d[seq_len(nc)],
    v = junk$v[, seq_len(nc)],
    feature.means = feature.means
  )

}
