#' Gene-specific Regularized Ordinary Least Squares fit statistics for supervised PCA
#'
#' @param x An \eqn{p \times n} predictor matrix.
#' @param y A response vector.
#' @param s0.perc Percentile of the standard error of the slope estimate to be
#'   used for regularization. The Default value of \code{NULL} will use the
#'   median of this distribution.
#'
#' @description Model statistics for Ordinary Least Squares (OLS) regression by
#'   gene.
#'
#' @return A list of OLS model statistics:
#' \itemize{
#'   \item{\code{tt} : }{The Student's \eqn{t} test statistic the slopes
#'      (\eqn{\beta}).}
#'   \item{\code{numer} : }{The estimate of \eqn{\beta}.}
#'   \item{\code{sd} : }{The standard error of the estimates for \eqn{\beta}
#'      (the standard error divided by the square root of Sxx).}
#'   \item{\code{fudge} : }{A regularization parameter. See Details for
#'      description.}
#' }
#'
#' @details This function calculates the Sxx, Syy, and Sxy sums from the gene-
#'   specific OLS models, then calculates estimates of the regression slopes for
#'   each gene and their corresponding regularized test statistics,
#'   \deqn{t = \hat{\beta} / (sd + e),}
#'   where \eqn{e} is a regularization parameter.
#'
#'   If \code{s0.perc} is \code{NULL}, then \eqn{e} is median of the \code{sd}
#'   values. Otherwise, \eqn{e} is set equal to \code{quantile(sd, s0.perc)}.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use SuperPCA_pVals() instead
#'   
#'   
#'   p <- 500
#'   n <- 50
#' 
#'   x_mat <- matrix(rnorm(n * p), nrow = p, ncol = n)
#'   time_int <- rpois(n, lambda = 365 * 2)
#'   
#'   olsTrain_fun(
#'     x = x_mat,
#'     y = time_int
#'   )
#'   
#'   
olsTrain_fun <- function(x, y, s0.perc = NULL){
  # browser()

  n <- length(y)
  xbar <- x %*% rep(1 / n, n)
  sxx <- ((x - as.vector(xbar)) ^ 2) %*% rep(1, n)
  sxy <- (x - as.vector(xbar)) %*% (y - mean(y))
  syy <- sum((y - mean(y)) ^ 2)

  estBetas <- sxy / sxx

  SSE <- syy - sxy ^ 2 / sxx
  std_Err <- sqrt(SSE / (n - 2))
  sd <- std_Err / sqrt(sxx)

  if(is.null(s0.perc)){
    fudge <- median(sd)
  } else {
    fudge <- ifelse(s0.perc >= 0, quantile(sd, s0.perc), 0)
  }

  tt <- estBetas / (sd + fudge)

  ###  Return  ###
  out <- list(tt = tt, numer = estBetas,
              sd = sd, fudge = fudge)
  out

}
