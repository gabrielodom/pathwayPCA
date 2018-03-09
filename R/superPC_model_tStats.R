#' Extract and Test Supervised PCs
#'
#' @description Identify significant features, extract PCs from those specific
#'    features to construct a data matrix, predict the response with this data
#'    matrix, and record the model fit statistic of this prediction.
#'
#' @param fit An object of class \code{superpc} returned by the function
#'    \code{superpc.train}
#' @param data A list of testing data:
#' \itemize{
#'   \item{x : }{A "tall" pathway data frame ($p_{path} * n$).}
#'   \item{y : }{A response vector corresponding to \code{type}}
#'   \item{censoring.status : }{If \code{type = "survival"}, the censoring
#'     indicator. Otherwise, \code{NULL}.}
#'   \item{featurenames : }{A character vector of the measured -omes in
#'     \code{x}.}
#'  }
#' @param n.threshold The number of bins into which to split the feature scores
#'    returned in the \code{fit} object.
#' @param threshold.ignore Calculate the model for feature scores above this
#'    percentile of the threshold. We have seen that the smalles threshold
#'    values (0\% - 40\%) largely have no effect on model t-scores. Defaults to
#'    0.00 (0\%).
#' @param n.PCs The number of PCs to extract from the significant pathway
#' @param min.features What is the smallest number of genes allowed in each
#'    pathway? This argument must be kept constant across all calls to this
#'    function which use the same pathway list. Defaults to 5
#' @param epsilon I'm not sure why this is important. It's called when comparing
#'    the absolute score values to each value of the threshold vector. Defaults
#'    to 10 ^ -6.
#'
#' @return A list containing:
#' \itemize{
#'   \item{thresholds : }{A labelled vector of quantile values of the score
#'     vector in the \code{fit} object.}
#'   \item{n.threshold : }{The number of splits to make in the score vector.}
#'   \item{scor : }{A matrix of model fit statistics. Each column is the
#'     threshold level of predictors allowed into the model, and each row is a
#'     PC included. Which genes are included in the matrix before PC extraction
#'     is governed by comparing their model score to the quantile value of the
#'     scores at each threshold value.}
#'   \item{tscor : }{A matrix of model t-statisics for each PC included (rows)
#'     at each threshold level (columns).}
#'   \item{type : }{Which model was called? Options are survival, regression,
#'     or binary.}
#' }
#'
#' @details See \url{https://web.stanford.edu/~hastie/Papers/spca_JASA.pdf}
#'   An issue, the number of thresholds at which to test (\code{n.threshold}),
#'   can be larger than the number of features to bin. This is why so many of
#'   the t-statistics are constant - because the model isn't changing.
#'
#' @export
#'
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats quantile
#' @importFrom survival coxph
#' @importFrom survival coxph.control
#' @importFrom survival Surv
#'
#' @examples
#'   NULL

superpc.st <- function(fit,
                       data,
                       n.threshold = 20,
                       threshold.ignore = 0.00,
                       n.PCs = 1,
                       min.features = 5,
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
    cur.tt1 <- cur.tt1[1:500]
    thresholds <- quantile(cur.tt1, probs = thresh_probs)

  }

  scor  <- array(NA, c(n.PCs, n.threshold))
  tscor <- array(NA, c(n.PCs, n.threshold))
  scor.preval <- matrix(NA, nrow = n.PCs, ncol = n.threshold)
  scor.lower  <- NULL
  scor.upper  <- NULL
  v.preval <- array(NA, c(n, n.PCs, n.threshold))

  kept.thresholds <- which(thresh_probs >= threshold.ignore)

  for(i in kept.thresholds){
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
    cur.v <- as.matrix(cur.v.all[, 1:n.PCs.eff])

    for(k in 1:n.PCs){

      switch(type,
             survival = {

               junk <- coxph(response ~ cur.v[, 1:k],
                             control = coxph.control(iter.max = 10))
               scor[k, i]  <- 2 * (junk$loglik[2] - junk$loglik[1])
               # Technically a z-score
               tscor[k, i] <- summary(junk)$coef[k, 4]

             },
             regression = {

               junk <- summary(lm(response ~ cur.v[, 1:k]))
               scor[k, i]  <- junk$fstat[1]
               tscor[k, i] <- junk$coef[k + 1, 3]

             },
             binary = {

               junk <- summary(glm(response ~ cur.v[, 1:k],
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


  junk <- list(thresholds = thresholds,
               n.threshold = n.threshold,
               scor = scor,
               tscor = tscor,
               type = type)
  class(junk) <- "superpc.st"
  return(junk)

}

