#' Gene-specific Generalized Linear Model fit statistics for supervised PCA
#'
#' @description Model statistics for Generalized Linear Model (GLM) regression
#'   by gene
#'
#' @param x An \eqn{p \times n} predictor matrix.
#' @param y A response vector.
#' @param family A description of the error distribution and link function to
#'   be used in the model. The default is \code{binomial(link = "logit")}.
#'
#' @return The slope coefficient from the GLM for each gene.
#'
#' @details While this function currently supports any GLM family from the
#'   \code{\link[stats]{family}} function, this function is only called in the
#'   model fitting step (via the internal \code{\link{superpc.train}}) function
#'   and not in the test statistic calculation step (in the
#'   \code{\link{superpc.st}} function). We would like to support Poisson
#'   regression through the \code{\link[stats]{glm}} function, as well as n-ary
#'   classification through \code{\link[nnet]{multinom}} and ordinal logistic
#'   regression through \code{\link[MASS]{polr}}.
#'
#' @keywords internal
#'
#' @importFrom stats binomial
#' @importFrom stats glm
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
#'   obs_logi <- sample(
#'     c(FALSE, TRUE),
#'     size = n,
#'     replace = TRUE,
#'     prob = c(0.2, 0.8)
#'   )
#' 
#'   glmTrain_fun(
#'     x = x_mat,
#'     y = obs_logi
#'   )
#'   
#'   
glmTrain_fun <- function(x, y, family = binomial){

  glmCoeffs <- function(x_mat){
    glm(y ~ x_mat, family = family)$coef[2]
  }

  tt <- apply(x, 1, glmCoeffs)


  return(list(tt = tt))
}
