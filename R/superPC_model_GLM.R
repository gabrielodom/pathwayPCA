#' Gene-Specific Generalized Linear Model Estimates
#'
#' @param x An \eqn{p \times n} predictor matrix
#' @param y A response vector
#' @param family a description of the error distribution and link function to
#'   be used in the model. The default is \code{binomial(link = "logit")}.
#'
#' @description Model statistics for Generalized Linear Model (GLM) regression
#'   by gene
#'
#' @return The slope coefficient from the GLM for each gene
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
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats poisson
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use superPCA_pVals() instead
glmTrain_fun <- function(x, y, family = binomial){

  glmCoeffs <- function(x_mat){
    glm(y ~ x_mat, family = family)$coef[2]
  }

  tt <- apply(x, 1, glmCoeffs)


  return(list(tt = tt))
}
