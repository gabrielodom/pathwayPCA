#' @title Parametric bootstrap and non-parametric permutations of a response vector or
#'    matrix
#'
#' @description Create a random parametric bootstrap sample or a permutation of
#'    the input response vector or matrix (for survival outcomes).
#'
#' @param response_vec The dependent vector to sample from. For survival
#'    response, this is the vector of event times. For regression or n-ary
#'    classification, this is the vector of responses.
#' @param event_vec The death / event observation indicator vector for survival
#'    response. This is coded as 0 for a right-censoring occurence and 1 for a
#'    recorded event.
#' @param respType What type of response has been supplied. Options are
#'    \code{"none"}, \code{"survival"}, \code{"regression"}, and
#'    \code{"categorical"}. Defaults to \code{"none"} to match the default
#'    \code{response = NULL} value.
#' @param parametric Should the random sample be taken using a parametric
#'    bootstrap sample? Defaults to \code{TRUE}.
#'
#' @return If \code{parametric = FALSE}, a permutation of the supplied response
#'    is returned (for AES-PCA). If \code{parametric = TRUE}, we return a
#'    parametric bootstrap sample of the response.
#'
#' @details The distributions (for \code{parametric = TRUE}) are Weibull for
#'    survival times, Normal for regression response, and n-ary Multinomial for
#'    categorical response. Distributional parameters are estimated with their
#'    maximum likelihood estimates. When \code{parametric = FALSE}, the response
#'    vector or survival matrix is randomly ordered by row. This option should
#'    only be used when called from the AES-PCA method.
#'
#' @keywords internal
#'
#' @importFrom survival Surv
#' @importFrom survival survreg
#' @importFrom stats rweibull
#' @importFrom stats runif
#'
#' @importFrom stats rnorm
#' @importFrom stats sd
#'
#' @importFrom methods as
#'
#'
#' @examples
#'   # DO NOT CALL THESE FUNCTIONS DIRECTLY.
#'   # Use AESPCA_pVals() or SuperPCA_pVals() instead
#'
#'   data("colon_pathwayCollection")
#'   data("colonSurv_df")
#'
#'   SampleResponses(
#'     response_vec = colonSurv_df$OS_time,
#'     event_vec = colonSurv_df$OS_event,
#'     respType = "survival"
#'   )
#'
#' @name RandomControlSample
#' @rdname permuteSamps
NULL



#' @export
#' @rdname permuteSamps
SampleResponses <- function(response_vec,
                            event_vec = NULL,
                            respType = c("survival",
                                         "regression",
                                         "categorical"),
                            parametric = TRUE){

  switch(respType,
         survival = {
           SampleSurv(
             response_vec, event_vec, parametric = parametric
           )
         },
         regression = {
           SampleReg(
             response_vec, parametric = parametric
           )
         },
         categorical = {
           SampleCateg(
             response_vec, parametric = parametric
           )
         })

}


#' @export
#' @rdname permuteSamps
SampleSurv <- function(response_vec,
                       event_vec,
                       parametric = TRUE){

  # browser()
  n <- length(event_vec)

  # Set up event time sample: parametric bootstrap or non-parametric permutation
  if(parametric){

    # Estimate the null model
    surv_obj <- Surv(response_vec, event_vec)
    null_mod <- survreg(surv_obj ~ 1, dist = "weibull")

    times_vec <- rweibull(
      n, shape = 1 / null_mod$scale, scale = exp(null_mod$coef)
    )

    # Randomly censor some of the observations
    pCensor <- mean(event_vec == 0)
    event_ind <- runif(n) >= pCensor
    for(m in 1:n){

      if(!event_ind[m]){
        times_vec[m] <- runif(1, min = 0, max = times_vec[m])
      }

    }

  } else {

    randIdx <- sample.int(n, n)
    times_vec <- response_vec[randIdx]
    event_ind <- event_vec[randIdx]

  }

  # Return the bootstrapped / permuted survival data
  list(response_vec = times_vec, event_vec = event_ind)

}



#' @export
#' @rdname permuteSamps
SampleReg <- function(response_vec,
                      parametric = TRUE){

  # browser()
  n <- length(response_vec)

  # Set up event time sample: parametric bootstrap or non-parametric permutation
  if(parametric){
    rand_vec <- rnorm(n, mean = mean(response_vec), sd = sd(response_vec))
  } else {
    rand_vec <- sample(response_vec)
  }

  # Return the bootstrapped / permuted regression response data
  rand_vec

}



#' @export
#' @rdname permuteSamps
SampleCateg <- function(response_vec,
                        parametric = TRUE){

  # browser()
  responseType <- class(response_vec)
  x <- as.factor(response_vec)

  # Set up event time sample: parametric bootstrap or non-parametric permutation
  if(parametric){

    n <- length(response_vec)
    classes <- unique(x)
    classProbs <- summary(x) / length(x)
    rand_vec <- sample(classes, size = n, replace = TRUE, prob = classProbs)

  } else {
    rand_vec <- sample(x)
  }

  # Return the bootstrapped / permuted regression response data in its original
  #   class
  if(responseType == "factor"){
    rand_vec
  } else {
    as(as.character(rand_vec), responseType)
  }

}
