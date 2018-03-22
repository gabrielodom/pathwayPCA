#' Calculate pathway-specific Student's \eqn{t}-scores from a null distribution
#'    for supervised PCA
#'
#' @description Randomly permute or parametrically resample the response vector
#'    before model analysis. Then extract principal components (PCs) from the
#'    gene pathway, and return the test statistics associated with the first
#'    \code{numPCs} principal components at a set of threshold values based on
#'    the permuted values of the response.
#'
#' @param pathway_vec A character vector of the measured -Omes in the chosen
#'    gene pathway. These should match a subset of the rownames of the gene
#'    array.
#' @param geneArray_df A "tall" pathway data frame (\eqn{p \times N}). Each
#'    subject or tissue sample is a column, and the rows are the -Ome
#'    measurements for that sample.
#' @param response_mat A response matrix corresponding to \code{responseType}.
#'    For \code{"regression"} and \code{"classification"}, this will be an
#'    \eqn{N \times 1} matrix of response values. For \code{"survival"}, this
#'    will be an \eqn{N \times 2} matrix with event times in the first column
#'    and observed event indicator in the second.
#' @param responseType A character string. Options are \code{"survival"},
#'    \code{"regression"}, and \code{"classification"}.
#' @param parametric Should the random sample be taken using a parametric
#'    bootstrap sample? Defaults to \code{FALSE}.
#' @param n.threshold The number of bins into which to split the feature scores
#'    in the \code{fit} object returned internally by the
#'    \code{\link{superpc.train}} function.
#' @param numPCs The number of PCs to extract from the pathway.
#' @param min.features What is the smallest number of genes allowed in each
#'    pathway? This argument must be kept constant across all calls to this
#'    function which use the same pathway list. Defaults to 3.
#'
#' @return A matrix with \code{numPCs} rows and \code{n.threshold} columns.
#'    The matrix values are model \eqn{t}-statisics for each PC included (rows)
#'    at each threshold level (columns).
#'
#' @details This is a wrapper function to call \code{\link{superpc.train}}
#'    and \code{\link{superpc.st}} after response sampling or permutation with
#'    the \code{\link{randomControlSample}} suite of functions. This response
#'    randomization will act as a null distribution against which to compare
#'    the results from the \code{\link{pathway_tScores}} function.
#'
#'    This wrapper is designed to facilitate apply calls (in parallel or
#'    serially) of these two functions over a list of gene pathways. When
#'    \code{numPCs} is equal to 1, we recommend using a simplify-style apply
#'    variant, such as \code{sapply} (shown in \code{\link[base]{lapply}}) or
#'    \code{parSapply} (shown in \code{\link[parallel]{clusterApply}}), then
#'    transposing the resulting matrix.
#'
#' @seealso \code{\link{pathway_tScores}}; \code{\link{randomControlSample}};
#'    \code{\link{superpc.train}}; \code{\link{superpc.st}}
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use superPCA_pVals() instead
pathway_tControl <- function(pathway_vec,
                             geneArray_df,
                             response_mat,
                             responseType = c("survival",
                                              "regression",
                                              "classification"),
                             parametric = FALSE,
                             n.threshold = 20,
                             numPCs = 1,
                             min.features = 3){
  # browser()

  data_ls <- switch(responseType,
    survival = {

      surv_ls <- sample_Survivalresp(response_vec = response_mat[, 1],
                                     censor_vec = response_mat[, 2],
                                     parametric = parametric)
      list(x = geneArray_df[pathway_vec, ],
           y = surv_ls$response_vec,
           censoring.status = surv_ls$censor_vec,
           featurenames = pathway_vec)

      },
    regression = {

      list(x = geneArray_df[pathway_vec, ],
           y = sample_Regresp(response_vec = response_mat,
                              parametric = parametric),
           featurenames = pathway_vec)

      },
    classification = {

      list(x = geneArray_df[pathway_vec, ],
           y = sample_Classifresp(response_vec = response_mat,
                                  parametric = parametric),
           featurenames = pathway_vec)

      }
    )

  train <- superpc.train(data_ls, type = responseType)

  st.obj <- superpc.st(fit = train,
                       data = data_ls,
                       n.PCs = numPCs,
                       min.features = min.features,
                       n.threshold = n.threshold)

  st.obj$tscor

}
