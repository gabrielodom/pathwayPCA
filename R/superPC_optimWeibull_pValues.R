#' Calculate the p-Values from an Optimal Mixture of Weibull Distributions
#'
#' @description Calculate the p-values of test statistics from a mixture of two
#'   Weibull Extreme Value distributions.
#'
#' @param tScore_vec A vector of the maximum absolute t-scores for each pathway
#'   when under the alternative model
#' @param pathwaySize_vec A vector of the number of genes in each pathway
#' @param optimParams_vec The NAMED vector of optimal mixture distribution
#'   parameters returned by the \code{\link{weibullMix_optimParams}} function.
#'
#' @return A named vector of the estimated raw p-values for each gene pathway
#'
#' @details The likelihood function is equation (4) Chen et al (2008): a mixture
#'  of two Gumbel Extreme Value pdfs, with mixing proportion p. The values mu_i
#'  and s_i, i = 1, 2, within the code of the function are placeholders for the
#'  mean and standard deviation, respectively.
#'
#'  @seealso \code{\link[stats]{optim}}
#'
#' @export
#'
#' @examples
#'    NULL

weibullMix_pValues <- function(tScore_vec,
                               pathwaySize_vec,
                               optimParams_vec){

  ###  Pathway Cardinality  ###
  # Still don't know what this does
  calc_anbn <- function(length_vec){

    logn <- log(length_vec)
    an1 <- sqrt(2 * logn)
    top <- log(4 * pi) + log(logn)
    bottom <- 2 * logn
    bn1 <- an1 * (1 - 0.5 * top / bottom)

    list(an = an1, bn = bn1)

  }
  abn_ls <- calc_anbn(pathwaySize_vec)


  ###  Density Values at Extremum  ###
  an_s1 <- abn_ls$an * optimParams_vec["s1"]
  arg_A <- -(tScore_vec - abn_ls$bn - optimParams_vec["mu1"]) * an_s1
  A <- 1 - exp(-exp(arg_A)) # in [0, 1]

  an_s2 <- abn_ls$an * optimParams_vec["s2"]
  arg_B <- (tScore_vec + abn_ls$bn + optimParams_vec["mu2"]) * an_s2
  B <- exp(-exp(arg_B)) # also in [0, 1]

  # These values will be between 0 and 1 if the "p" value is in [0, 1]
  tt1 <- optimParams_vec["p"] * A + (1 - optimParams_vec["p"]) * B
  tt2 <- 1 - tt1


  ###  Return p-Values  ###
  apply(cbind(tt1, tt2), MARGIN = 1, FUN = min)

}
