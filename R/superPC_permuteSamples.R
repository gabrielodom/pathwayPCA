#' Parametric Bootstrap and Non-parametric Permutations of a Response
#'
#' Create a random parametric bootstrap sample or a permutation of
#'    the input response vector or matrix (for survival outcomes).
#'
#' @param response_vec The dependent vector to sample from. For survival
#'    response, this is the vector of event times. For regression or n-ary
#'    classification, this is the vector of responses.
#' @param censor_vec The censoring indicator vector for survival response. This
#'    is coded as 1 for a right-censoring occurence and 0 for a recorded event.
#' @param parametric Should the random sample be taken using a parametric
#'    bootstrap sample? Defaults to \code{FALSE}.
#'
#' @details The distributions (for \code{parametric = TRUE}) are Weibull for
#'    survival times, Normal for regression, and n-ary Multinomial for
#'    classification. Distributional parameters are estimated with their maximum
#'    likelihood estimates. When \code{parametric = FALSE}, the response vector
#'    or survival matrix is simply permuted by row.
#'
#' @name randomControlSample
NULL



#'
#' @importFrom survival Surv
#' @importFrom survival survreg
#' @importFrom stats rweibull
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#'    NULL
#' @rdname randomControlSample
sample_Survivalresp <- function(response_vec,
                                censor_vec,
                                parametric = FALSE){

  # browser()
  n <- length(censor_vec)

  # Set up event time sample: parametric bootstrap or non-parametric permutation
  if(parametric){

    # Estimate the null model
    surv_obj <- survival::Surv(response_vec, censor_vec)
    null_mod <- survival::survreg(surv_obj ~ 1, dist = "weibull")

    times_vec <- rweibull(n,
                          shape = 1 / null_mod$scale,
                          scale = exp(null_mod$coef))

    # Randomly censor some of the observations
    pCensor <- mean(censor_vec == 1)
    censor_ind <- runif(n) < pCensor
    for(m in 1:n){

      if(censor_ind[m]){
        times_vec[m] <- runif(1, min = 0, max = times_vec[m])
      }

    }

  } else {

    randIdx <- sample.int(n, n)
    times_vec <- response_vec[randIdx]
    censor_ind <- censor_vec[randIdx]

  }

  # Return the bootstrapped / permuted survival data
  list(response_vec = times_vec, censor_vec = censor_ind)

}




#'
#' @importFrom stats rnorm
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#'    NULL
#' @rdname randomControlSample
sample_Regresp <- function(response_vec,
                           parametric = FALSE){

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



#'
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#'    NULL
#' @rdname randomControlSample
sample_Classifresp <- function(response_vec,
                               parametric = FALSE){

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
