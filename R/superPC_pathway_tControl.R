#' Calculate pathway-specific Student's \eqn{t}-scores from a null distribution
#'    for supervised PCA
#'
#' @description Parametrically resample the response vector before model
#'    analysis. Then extract principal components (PCs) from the gene pathway,
#'    and return the test statistics associated with the first \code{numPCs}
#'    principal components at a set of threshold values based on the resampled
#'    values of the response.
#'
#' @param pathway_vec A character vector of the measured -Omes in the chosen
#'    gene pathway. These should match a subset of the rownames of the gene
#'    array.
#' @param geneArray_df A "tall" pathway data frame (\eqn{p \times N}). Each
#'    subject or tissue sample is a column, and the rows are the -Ome
#'    measurements for that sample.
#' @param response_mat A response matrix corresponding to \code{responseType}.
#'    For \code{"regression"} and \code{"categorical"}, this will be an
#'    \eqn{N \times 1} matrix of response values. For \code{"survival"}, this
#'    will be an \eqn{N \times 2} matrix with event times in the first column
#'    and observed event indicator in the second.
#' @param responseType A character string. Options are \code{"survival"},
#'    \code{"regression"}, and \code{"categorical"}.
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
#'    and \code{\link{superpc.st}} after response parametric bootstrapping with
#'    the \code{\link{RandomControlSample}} suite of functions. This response
#'    sampling will act as a null distribution against which to compare
#'    the results from the \code{\link{pathway_tScores}} function.
#'
#'    This wrapper is designed to facilitate apply calls (in parallel or
#'    serially) of these two functions over a list of gene pathways. When
#'    \code{numPCs} is equal to 1, we recommend using a simplify-style apply
#'    variant, such as \code{sapply} (shown in \code{\link[base]{lapply}}) or
#'    \code{parSapply} (shown in \code{\link[parallel]{clusterApply}}), then
#'    transposing the resulting matrix.
#'
#' @seealso \code{\link{pathway_tScores}}; \code{\link{RandomControlSample}};
#'    \code{\link{superpc.train}}; \code{\link{superpc.st}}
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use SuperPCA_pVals() instead
#'
#'   data("colon_pathwayCollection")
#'   data("colonSurv_df")
#'
#'   colon_OmicsSurv <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "surv"
#'   )
#'
#'   asthmaGenes_char <-
#'     getTrimPathwayCollection(colon_OmicsSurv)[["KEGG_ASTHMA"]]$IDs
#'   resp_mat <- matrix(
#'     c(getEventTime(colon_OmicsSurv), getEvent(colon_OmicsSurv)),
#'     ncol = 2
#'   )
#'
#'   pathway_tControl(
#'     pathway_vec = asthmaGenes_char,
#'     geneArray_df = t(getAssay(colon_OmicsSurv)),
#'     response_mat = resp_mat,
#'     responseType = "survival"
#'   )
#'
pathway_tControl <- function(pathway_vec,
                             geneArray_df,
                             response_mat,
                             responseType = c("survival",
                                              "regression",
                                              "categorical"),
                             n.threshold = 20,
                             numPCs = 1,
                             min.features = 3){
  # browser()

  sampResp <- SampleResponses(
    response_vec = response_mat[, 1],
    event_vec = response_mat[, 2],
    respType = responseType,
    parametric = TRUE
  )

  data_ls <- switch(responseType,
    survival = {
      list(
        x = geneArray_df[pathway_vec, ],
        y = sampResp$response_vec,
        censoring.status = sampResp$event_vec,
        featurenames = pathway_vec
      )
    },
    regression = {
      list(
        x = geneArray_df[pathway_vec, ],
        y = sampResp,
        featurenames = pathway_vec
      )
    },
    categorical = {
      list(
        x = geneArray_df[pathway_vec, ],
        y = sampResp,
        featurenames = pathway_vec
      )
    })

  train <- superpc.train(data_ls, type = responseType)

  st.obj <- superpc.st(
    fit = train,
    data = data_ls,
    n.PCs = numPCs,
    min.features = min.features,
    n.threshold = n.threshold
  )

  st.obj$tscor

}
