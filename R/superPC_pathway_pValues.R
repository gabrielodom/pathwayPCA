#' Calculate the p-Values from a Mixture if Weibull Extreme Value Distributions
#'
#' @description Calculate the pathway-specific p-values and their associated
#'    False Discovery Rates (FDR).
#'
#' @param optimParams_vec A named vector of the estimated values for the
#'    parameters which minimize the likelihood as returned by the function
#'    \code{\link{weibullMix_optimParams}}.
#' @param max_tScores_vec A vector of the maximum absolute t-scores for each
#'    pathway when under the alternative model (the response vector as is).
#' @param genelist_ls A list of three elements:
#'   \itemize{
#'     \item{pathways : }{A list of character vectors such that each vector
#'        contains the ID numbers (as a character) of the individual genes
#'        within that pathway as a vector of character strings.}
#'     \item{TERMS : }{A character vector containing the names of the gene
#'        pathways.}
#'     \item{setsize : }{A named integer vector containing the number of genes
#'        in each gene pathway.}
#'   }
#' @param FDRadjust Should the p-values be adjusted for multiple comparisons?
#'    Defaults to \code{TRUE}.
#' @param multTestProc If the p-values should be adjusted, which procedure will
#'    you use? Options are passed to the \code{\link{adjustRaw_pVals}} function.
#'    Specify multiple procedures via \code{c(...)}. Defaults to "BH".
#'
#' @return A data frame
#'
#' @details This function takes in the optimal parameters returned by the
#'    \code{\link{weibullMix_optimParams}} function, the maximum t-Scores for
#'    each gene pathway, and the list of gene pathway information. This function
#'    will calculate the p-value for each t-Score given the Gumbel Extreme Value
#'    mixture distribution parametrized by the values of \code{optimParams}. If
#'    requested, this function will also calculate the FDR associated with all
#'    pathway p-values via the Benjamini & Hochberg (1995) step-up FDR-
#'    controlling procedure as implemented in the \code{multtest} package.
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use superPCA_pVals() instead
pathway_pValues <- function(optimParams_vec,
                            max_tScores_vec,
                            genelist_ls,
                            FDRadjust = TRUE,
                            multTestProc = "BH"){
  # browser()

  ###  Pathway Cardinality  ###
  # Don't know what this does, but it's also in the weibullMix_optimParams
  #   function
  calc_anbn <- function(length_vec){

    logn <- log(length_vec)
    an1 <- sqrt(2 * logn)
    top <- log(4 * pi) + log(logn)
    bottom <- 2 * logn
    bn1 <- an1 * (1 - 0.5 * top / bottom)

    list(an = an1, bn = bn1)

  }
  abn_ls <- calc_anbn(genelist_ls$setsize)


  ###  Weibull Mixture p-Value Calculation  ###
  newP_fun <- function(tScore_vec, optimParams_vec, abCounts_ls){

    an_s1 <- abCounts_ls$an * optimParams_vec["s1"]
    arg_A <- -(tScore_vec - abCounts_ls$bn - optimParams_vec["mu1"]) * an_s1
    A <- 1 - exp(-exp(arg_A)) # in [0, 1]

    an_s2 <- abCounts_ls$an * optimParams_vec["s2"]
    arg_B <- (tScore_vec + abCounts_ls$bn + optimParams_vec["mu2"]) * an_s2
    B <- exp(-exp(arg_B)) # also in [0, 1]

    # These values will be between 0 and 1 if the "p" value is in [0, 1]
    tt1 <- optimParams_vec["p"] * A + (1 - optimParams_vec["p"]) * B
    tt2 <- 1 - tt1

    apply(cbind(tt1, tt2), MARGIN = 1, FUN = min)

  }

  newp <- newP_fun(tScore_vec = max_tScores_vec,
                   optimParams_vec = optimParams_vec,
                   abCounts_ls = abn_ls)
  ntest <- data.frame(goterms = names(genelist_ls$pathways),
                      setsize = genelist_ls$setsize,
                      rawp = newp)
  rownames(ntest) <- NULL


  ###  p-Value Adjustment  ###

  if(FDRadjust){
    # browser()

    adjustedP <- adjustRaw_pVals(ntest$rawp, multTestProc)

    ntest <- cbind(ntest, adjustedP[, -1, drop = FALSE])
    ntest$terms <- genelist_ls$TERMS

  }


  ###  Return  ###
  ntest

}
