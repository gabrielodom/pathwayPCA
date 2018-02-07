#' Calculate Pathway-Specific Student's t Scores
#'
#' @param pathway_vec A character vector of the measured -omes in the chosen
#'    gene pathway. These should match a subset of the rownames of the gene
#'    array.
#' @param geneArray_df A "tall" pathway data frame ($p * n$). Each subject or
#'    tissue sample is a column, and the rows are the -ome measurements for that
#'    sample.
#' @param response_mat A response matrix corresponding to \code{responseType}.
#'    For \code{"regression"} and \code{"classification"}, this will be an
#'    $n * 1$ matrix of response values. For \code{"survival"}, this will be an
#'    $n * 2$ matrix with event times in the first column and censoring
#'    indicator in the second.
#' @param responseType A character string. Options are \code{"survival"},
#'    \code{"regression"}, and \code{"classification"}.
#' @param n.threshold The number of bins into which to split the feature scores
#'    in the \code{fit} object returned internally by the \code{superpc.train}
#'    function.
#' @param numPCs The number of PCs to extract from the significant pathway.
#' @param min.features What is the smallest number of genes allowed in each
#'    pathway? This argument must be kept constant across all calls to this
#'    function which use the same pathway list. Defaults to 5.
#'
#' @description Extract principal components from the gene pathway, and return
#'    the test statistics associated with the first \code{numPCs} principal
#'    components at a set of threshold values.
#'
#' @return A matrix with \code{numPCs} rows and \code{n.threshold} columns:
#'    model t-statisics for each PC included (rows) at each threshold level
#'    (columns).
#'
#' @details This is a wrapper function to call \code{\link{superpc.train}}
#'    and \code{\link{superpc.st}}. This wrapper is designed to facilitate apply
#'    calls (in parallel or serially) of these two functions over a list of gene
#'    pathways. When \code{numPCs} is 1, we recommend using a simplify-style
#'    apply variant, such as \code{sapply} (\code{\link[base]{lapply}}) or
#'    \code{parSapply} (\code{\link[parallel]{clusterApply}}), then transposing
#'    the resulting matrix.
#'
#' @seealso \code{\link{superpc.train}}; \code{\link{superpc.st}}
#'
#' @export
#'
#' @examples
#'    NULL
pathway_tScores <- function(pathway_vec,
                            geneArray_df,
                            response_mat,
                            responseType = c("survival",
                                             "regression",
                                             "classification"),
                            n.threshold = 20,
                            numPCs = 1,
                            min.features = 5){
  # browser()

  data_ls <- switch(responseType,
                    survival = {
                      list(x = geneArray_df[pathway_vec, ],
                           y = response_mat[, 1],
                           censoring.status = response_mat[, 2],
                           featurenames = pathway_vec)
                    },
                    regression = {
                      list(x = geneArray_df[pathway_vec, ],
                           y = response_mat,
                           featurenames = pathway_vec)
                    },
                    classification = {
                      list(x = geneArray_df[pathway_vec, ],
                           y = response_mat,
                           featurenames = pathway_vec)
                    })

  train <- superpc.train(data_ls, type = responseType)

  st.obj <- superpc.st(fit = train,
                       data = data_ls,
                       n.PCs = numPCs,
                       min.features = min.features,
                       n.threshold = n.threshold)

  st.obj$tscor

}
