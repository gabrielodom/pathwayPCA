#' Calculate the Optimal Parameters for a Mixture of Weibull Distributions
#'
#' @description Calculate the parameters which minimise the likelihood from a
#'    mixture of two Weibull Extreme Value distributions.
#'
#' @param max_tControl_vec A vector of the maximum absolute t-scores for each
#'    pathway when under the null model (the response vector has been randomly
#'    assigned or parametrically bootstrapped).
#' @param pathwaySize_vec A vector of the number of genes in each pathway
#' @param initialVals A vector of initial values for the Weibull parameters.
#'    The values are
#'    \itemize{
#'      \item{p : }{The mixing proportion between the Gumbel minimum and Gumbel
#'         maximum distributions. Defaults to 0.5.}
#'      \item{mu1 : }{The mean of the first distribution. Defaults to 1.}
#'      \item{s1 : }{The standard deviation (precision?) of the first
#'         distribution. Defaults to 0.5.}
#'      \item{mu2 : }{The mean of the second distribution. Defaults to 1.}
#'      \item{s2 : }{The standard deviation (precision?) of the second
#'         distribution. Defaults to 0.5.}
#'    }
#' @param optimMethod Which numerical optimization routine to pass to the
#'    \code{optim} function. Defaults to "L-BFGS-B", which allows for lower and
#'    upper bound constraints. When this option is selected, lower and upper
#'    bounds for ALL parameters must be supplied.
#' @param lowerBD A vector of the lower bounds on the \code{initialVals}.
#' @param upperBD A vector of the upper bounds on the \code{initialVals}.
#'
#' @return A named vector of the estimated values for the parameters which
#'    minimize the likelihood.
#'
#' @details The likelihood function is equation (4) Chen et al (2008): a mixture
#'  of two Gumbel Extreme Value pdfs, with mixing proportion p. The values mu_i
#'  and s_i, i = 1, 2, within the code of the function are placeholders for the
#'  mean and standard deviation, respectively.
#'
#'  A computational note: the \code{"L-BFGS-B"} option within the
#'  \code{\link[stats]{optim}} function requires a bounded function or
#'  likelihood. We therefore replaced \code{Inf} with $10 ^ 200$ in the
#'  check for boundedness. As we are attempting to minimise the likelihood, this
#'  maximum machine value is effectively Infinity.
#'
#'  See \url{https://doi.org/10.1093/bioinformatics/btn458} for more information
#'
#'  @seealso \code{\link[stats]{optim}}
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats deriv
#' @importFrom stats optim
#'
#' @examples
#'    NULL

weibullMix_optimParams <- function(max_tControl_vec,
                                   pathwaySize_vec,
                                   initialVals = c(p = 0.5,
                                                   mu1 = 1, s1 = 0.5,
                                                   mu2 = 1, s2 = 0.5),
                                   optimMethod = "L-BFGS-B",
                                   lowerBD = c(0, -Inf, 0, -Inf, 0),
                                   upperBD = c(1, Inf, Inf, Inf, Inf)){

  ###  Pathway Cardinality  ###
  # Don't know what this does
  calc_anbn <- function(length_vec){

    logn <- log(length_vec)
    an1 <- sqrt(2 * logn)
    top <- log(4 * pi) + log(logn)
    bottom <- 2 * logn
    bn1 <- an1 * (1 - 0.5 * top / bottom)

    list(an = an1, bn = bn1)

  }
  abn_ls <- calc_anbn(pathwaySize_vec)


  ###  The Likelihood Function  ###
  # The formulation here is different from the paper in the following:
  #   1. They define zi = (t - mu_i) / sigma_i
  #   2. They divide p and (1 - p) by sigma 1 and 2, respectively.
  # The values for s are not used appropriately: the s values are used more as
  #   precision (being multiplied instead of used as a divisor).
  gumbelMixture <- function(p_vec, maxt_vec, an_vec, bn_vec){

    z1 <- (maxt_vec - bn_vec - p_vec["mu1"]) * an_vec * p_vec["s1"]
    z2 <- (maxt_vec + bn_vec + p_vec["mu2"]) * an_vec * p_vec["s2"]
    e  <- (p_vec["p"] * an_vec * p_vec["s1"]) *
      exp(-z1 - exp(-z1)) +
      ((1 - p_vec["p"]) * an_vec * p_vec["s2"]) *
      exp(z2 - exp(z2))

    ifelse(test = any(e <= 0), yes = 10 ^ 200, no = -sum(log(e)))

  }


  ###  The Score Function  ###
  gumbelMix_score <- function(p_vec, maxt_vec, an_vec, bn_vec){

    innard_formula <- as.formula(~ -log(
      (p * an * s1) *
        exp(-((x - bn - mu1) * an * s1) - exp(-(x - bn - mu1) * an * s1)) +
        ((1 - p) * an * s2) *
        exp(((x + bn + mu2) * an * s2) - exp((x + bn + mu2) * an * s2))
    )
    )

    lmix_fun <- deriv(innard_formula,
                      c("p", "mu1", "s1", "mu2", "s2"),
                      function(x, an, bn, p, mu1, s1, mu2, s2) NULL)

    gradient <- lmix_fun(x = maxt_vec,
                         an = an_vec, bn = bn_vec,
                         p = p_vec["p"],
                         mu1 = p_vec["mu1"], s1 = p_vec["s1"],
                         mu2 = p_vec["mu2"], s2 = p_vec["s2"])
    # Extract the gradient matrix from the gradient vector, then find the
    #   column sums of that matrix
    colSums(attr(gradient,"gradient"))

  }


  ###  Constrained Optimization of the Likelihood  ###
  # Put a constraint on the bounds of p to avoid a global minimum where p is
  #   not an element of [0,1]. Otherwise, the solution will be p ~= -8k.
  pOptim <- optim(par = initialVals,
                  fn = gumbelMixture,
                  gr = gumbelMix_score,
                  maxt_vec = max_tControl_vec,
                  an_vec = abn_ls$an,
                  bn_vec = abn_ls$bn,
                  method = optimMethod,
                  lower = lowerBD,
                  upper = upperBD)$par

  pOptim

}
