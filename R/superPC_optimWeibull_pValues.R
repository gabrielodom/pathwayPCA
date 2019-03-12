#' Calculate the \eqn{p}-values from an optimal mixture of Weibull Extreme Value
#'    distributions for supervised PCA
#'
#' @description Calculate the \eqn{p}-values of test statistics from a mixture
#'    of two Weibull Extreme Value distributions.
#'
#' @param tScore_vec A vector of the maximum absolute \eqn{t}-scores for each
#'    pathway (returned by the \code{\link{pathway_tScores}} function) when
#'    under the alternative model.
#' @param pathwaySize_vec A vector of the number of genes in each pathway.
#' @param optimParams_vec The \emph{NAMED} vector of optimal Weibull Extreme
#'    Value mixture distribution parameters returned by the
#'    \code{\link{OptimGumbelMixParams}} function.
#'
#' @return A named vector of the estimated raw \eqn{p}-values for each gene
#'    pathway.
#'
#' @details The likelihood function is equation (4) in Chen et al (2008): a
#'    mixture of two Gumbel Extreme Value probability density functions, with
#'    mixing proportion \eqn{p}. Within the code of this function, the values
#'    \code{mu1}, \code{mu2} and \code{s1}, \code{s2} are placeholders for the
#'    mean and precision, respectively.
#'
#'    See \url{https://doi.org/10.1093/bioinformatics/btn458} for more
#'    information.
#'
#' @seealso \code{\link{OptimGumbelMixParams}}; \code{\link{pathway_tScores}};
#'    \code{\link{SuperPCA_pVals}}
#'
#' @keywords internal
#'
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use SuperPCA_pVals() instead.
#'
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colon_pathwayCollection")
#'   n_int <- lengths(colon_pathwayCollection$pathways)
#'
#'
#'   ###  Simulate Maximum Absolute Control t-Values  ###
#'   # The SuperPCA algorithm defaults to 20 threshold values; the example
#'   #   pathway collection has 15 pathways.
#'   t_mat <- matrix(rt(15 * 20, df = 5), nrow = 15)
#'
#'   absMax <- function(vec){
#'     vec[which.max(abs(vec))]
#'   }
#'   tAbsMax_num <- apply(t_mat, 1, absMax)
#'
#'
#'   ###  Calculate Optimal Parameters for the Gumbel Distribution  ###
#'   optParams_num <- OptimGumbelMixParams(
#'     max_tControl_vec = tAbsMax_num,
#'     pathwaySize_vec = n_int
#'   )
#'
#'
#'   ###  Simulate Maximum Absolute t-Values  ###
#'   tObs_mat <- matrix(rt(15 * 20, df = 3), nrow = 15)
#'   tObsAbsMax_num <- apply(tObs_mat, 1, absMax)
#'
#'
#'   ###  Calculate Observed-t-score p-Values  ###
#'   GumbelMixpValues(
#'     tScore_vec = tObsAbsMax_num,
#'     pathwaySize_vec = n_int,
#'     optimParams_vec = optParams_num
#'   )
#' }
#'
GumbelMixpValues <- function(tScore_vec,
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
  t_1 <- (tScore_vec - abn_ls$bn - optimParams_vec["mu1"]) * an_s1
  A <- exp(-exp(-t_1)) # in [0, 1]

  an_s2 <- abn_ls$an * optimParams_vec["s2"]
  t_2 <- (tScore_vec + abn_ls$bn + optimParams_vec["mu2"]) * an_s2
  B <- 1 - exp(-exp(t_2)) # also in [0, 1]

  # These values will be between 0 and 1 if the "p" value is in [0, 1]
  tt1 <- optimParams_vec["p"] * A + (1 - optimParams_vec["p"]) * B
  tt2 <- 1 - tt1


  ###  Return p-Values  ###
  apply(cbind(tt1, tt2), MARGIN = 1, FUN = min)

}
