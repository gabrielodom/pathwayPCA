#' Extract and test principal components from supervised PCA
#'
#' @description Identify \eqn{p_{path}} significant features, extract principal
#'    components (PCs) from those specific features to construct a data matrix,
#'    predict the response with this data matrix, and record the model fit
#'    statistic of this prediction.
#'
#' @param fit An object of class \code{superpc} returned by the function
#'    \code{\link{superpc.train}}.
#' @param data A list of test data:
#' \itemize{
#'   \item{\code{x} : }{A "tall" pathway data frame (\eqn{p_{path} \times N}).}
#'   \item{\code{y} : }{A response vector corresponding to \code{type}.}
#'   \item{\code{censoring.status} : }{If \code{type = "survival"}, the
#'      censoring indicator (\eqn{1 - } the observed event indicator).
#'      Otherwise, \code{NULL}.}
#'   \item{\code{featurenames} : }{A character vector of the measured -Omes in
#'     \code{x}.}
#'  }
#' @param n.threshold The number of bins into which to split the feature scores
#'    returned in the \code{fit} object.
#' @param threshold.ignore Calculate the model for feature scores above this
#'    percentile of the threshold. We have observed that the smallest threshold
#'    values (0\% - 40\%) largely have no effect on model \eqn{t}-scores.
#'    Defaults to 0.00 (0\%).
#' @param n.PCs The number of PCs to extract from the pathway.
#' @param min.features What is the smallest number of genes allowed in each
#'    pathway? This argument must be kept constant across all calls to this
#'    function which use the same pathway list. Defaults to 3.
#' @param epsilon I'm not sure why this is important. It's called when comparing
#'    the absolute score values to each value of the threshold vector. Defaults
#'    to \eqn{10^{-6}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{thresholds} : }{A labelled vector of quantile values of the
#'      score vector in the \code{fit} object.}
#'   \item{\code{n.threshold} : }{The number of splits to make in the score
#'      vector.}
#'   \item{\code{scor} : }{A matrix of model fit statistics. Each column is the
#'     threshold level of predictors allowed into the model, and each row is a
#'     PC included. Which genes are included in the matrix before PC extraction
#'     is governed by comparing their model score to the quantile value of the
#'     scores at each threshold value.}
#'   \item{\code{tscor} : }{A matrix of model \eqn{t}-statisics for each PC
#'      included (rows) at each threshold level (columns).}
#'   \item{\code{type} : }{Which model was called? Options are survival,
#'      regression, or binary.}
#' }
#'
#' @details NOTE: the number of thresholds at which to test (\code{n.threshold})
#'   can be larger than the number of features to bin. This will result in
#'   constant \eqn{t}-statistics for the first few bins because the model isn't
#'   changing.
#'
#'   See \url{https://web.stanford.edu/~hastie/Papers/spca_JASA.pdf}.
#'
#' @seealso \code{\link{superpc.train}}; \code{\link{SuperPCA_pVals}}
#'
#' @keywords internal
#'
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats quantile
#' @importFrom survival coxph
#' @importFrom survival coxph.control
#' @importFrom survival Surv
#'
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use SuperPCA_pVals() instead
#'   
#' \dontrun{
#'   data("colon_pathwayCollection")
#'   data("colonSurv_df")
#'   
#'   colon_OmicsSurv <- CreateOmics(
#'     assayData_df = colonSurv_df[,-(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "surv"
#'   )
#'   
#'   asthmaGenes_char <-
#'     getTrimPathwayCollection(colon_OmicsSurv)[["KEGG_ASTHMA"]]$IDs
#'     
#'   data_ls <- list(
#'     x = t(getAssay(colon_OmicsSurv))[asthmaGenes_char, ],
#'     y = getEventTime(colon_OmicsSurv),
#'     censoring.status = getEvent(colon_OmicsSurv),
#'     featurenames = asthmaGenes_char
#'   )
#'   
#'   superpcFit <- superpc.train(
#'     data = data_ls,
#'     type = "surv"
#'   )
#'   
#'   superpc.st(
#'     fit = superpcFit,
#'     data = data_ls
#'   )
#' } 
#'   
superpc.st <- function(fit,
                       data,
                       n.threshold = 20,
                       threshold.ignore = 0.00,
                       n.PCs = 1,
                       min.features = 3,
                       epsilon = 0.000001){
  # browser()

  type <- fit$type
  n <- ncol(data$x)
  p <- nrow(data$x)

  # This code had the survival object defined n.threshold * n.PCs
  #   times. This is unnecessary. Extract that code to here
  if(type == "survival"){
    response <- Surv(data$y, data$censoring.status)
  } else {
    response <- data$y
  }

  if(n.PCs > min.features){
    cat("Max # of components is min.features", fill = TRUE)
  }
  n.PCs <- min(min.features, n.PCs)

  cur.tt <- fit$feature.scores
  upper <- quantile(abs(cur.tt), 1 - (min.features / p))
  thresh_probs <- 0:(n.threshold - 1) / n.threshold

  if(length(cur.tt) < 500){

    if(quantile(abs(cur.tt), (n.threshold - 1) / n.threshold) > upper){

      cur.tt1 <- cur.tt
      cur.tt1[cur.tt1 > upper] <- upper

    } else {
      cur.tt1 <- cur.tt
    }

    thresholds <- quantile(cur.tt1, probs = thresh_probs)

  } else {

    cur.tt1 <- sort(cur.tt, decreasing = TRUE)
    cur.tt1 <- cur.tt1[seq_len(500)]
    thresholds <- quantile(cur.tt1, probs = thresh_probs)

  }

  scor  <- array(NA, c(n.PCs, n.threshold))
  tscor <- array(NA, c(n.PCs, n.threshold))
  scor.preval <- matrix(NA, nrow = n.PCs, ncol = n.threshold)
  scor.lower  <- NULL
  scor.upper  <- NULL
  v.preval <- array(NA, c(n, n.PCs, n.threshold))

  kept.thresholds <- which(thresh_probs >= threshold.ignore)

  # browser()
  for(i in kept.thresholds){
    # For many pathways, the first 10-15 of these thresholds yield identical
    #   results for the PCs 1 and 2. Additionally, all thresholds but the last
    #   two yield rather similar results. The last two are nearly spherical:
    #   the eigenvalues are within 30% of each other. The first 18 thresholds
    #   have the first eigenvalue at twice the size of the second
    # browser()

    # cat(i)
    cur.features <- (abs(cur.tt) + epsilon > thresholds[i])
    cur.svd <- mysvd(data$x[cur.features, ], n.components = n.PCs)
    xtemp <- data$x[cur.features, , drop = FALSE]
    xtemp <- t(scale(t(xtemp),
                     center = cur.svd$feature.means,
                     scale = FALSE))
    cur.v.all <- scale(t(xtemp) %*% cur.svd$u,
                       center = FALSE,
                       scale = cur.svd$d)
    n.PCs.eff <- min(sum(cur.features), n.PCs)
    cur.v <- as.matrix(cur.v.all[, seq_len(n.PCs.eff)])

    for(k in seq_len(n.PCs)){

      switch(type,
             survival = {

               junk <- coxph(response ~ cur.v[, seq_len(k)],
                             control = coxph.control(iter.max = 10))
               scor[k, i]  <- 2 * (junk$loglik[2] - junk$loglik[1])
               # Technically a z-score
               tscor[k, i] <- summary(junk)$coef[k, 4]

             },
             regression = {

               junk <- summary(lm(response ~ cur.v[, seq_len(k)]))
               scor[k, i]  <- junk$fstat[1]
               tscor[k, i] <- junk$coef[k + 1, 3]

             },
             binary = {

               junk <- summary(glm(response ~ cur.v[, seq_len(k)],
                                   family = binomial(link = "logit")))
               scor[k, i]  <- junk$null.deviance - junk$deviance
               tscor[k, i] <- junk$coef[k + 1, 3]

             },
             n_ary = {
               stop("Multinomial Regression not currently implemented.")
             },
             ordered = {
               stop("Ordered Logistic Regression not currently implemented.")
             })

    } # END for k

  } # END for i


  ###  Return the SuperPCA Results for the "Best" t-Score  ###
  # Find the most t-statistic(s), then pick the one that corresponds to the most
  #   extreme threshold
  bestT_idx <- which(abs(tscor[1,]) == max(abs(tscor[1,])))
  bestFeatures <- abs(cur.tt) + epsilon > thresholds[max(bestT_idx)]
  bestSVD <- mysvd(data$x[bestFeatures, ], n.components = n.PCs)
  bestLoadings_mat <- diag(bestSVD$d, ncol = length(bestSVD$d)) %*% t(bestSVD$u)
  rownames(bestLoadings_mat) <- paste0("PC", seq_len(n.PCs))
  colnames(bestLoadings_mat) <- names(bestSVD$feature.means)
  bestPCs_mat <- as.matrix(bestSVD$v, ncol = n.PCs)


  ###  Return  ###
  out_ls <- list(
    thresholds = thresholds,
    n.threshold = n.threshold,
    scor = scor,
    tscor = tscor,
    PCs_mat = bestPCs_mat,
    Loadings_mat = bestLoadings_mat,
    type = type
  )
  class(out_ls) <- "superpc.st"
  out_ls

}

