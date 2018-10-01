#' Train a supervised PCA model
#'
#' @description Computes feature scores for \eqn{p_{path}} features of a pathway
#'    via a linear model fit.
#'
#' @param data A list of test data:
#' \itemize{
#'   \item{\code{x} : }{A "tall" pathway data frame (\eqn{p_{path} \times N}).}
#'   \item{\code{y} : }{A response vector corresponding to \code{type}.}
#'   \item{\code{censoring.status} : }{If \code{type = "survival"}, the
#'      censoring indicator (\eqn{1 - } the observed event indicator. Otherwise,
#'      \code{NULL}.}
#'   \item{\code{featurenames} : }{A character vector of the measured -Omes in
#'     \code{x}.}
#'  }
#' @param type What model relates \code{y} and \code{x}? Options are
#'    \code{"survival"}, \code{"regression"}, or \code{"categorical"}.
#' @param s0.perc A stabilization parameter on the interval \eqn{[0,1]}. This is
#'    an internal argument to each of the called functions. The default value is
#'    \code{NULL} to ensure an appropriate value is determined internally.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{feature.scores} : }{The scaled \eqn{p}-dimensional score
#'      vector: each value has been divided by its respective standard deviation
#'      plus epsilon (governed by \code{s0.perc}). \code{NA} values returned by
#'      the logistic model are replaced with 0.}
#'   \item{\code{type} : }{The argument for \code{type}.}
#'   \item{\code{s0.perc} : }{The user-supplied value of \code{s0.perc}, or the
#'     internally-calculated default value from the chosen model.}
#'   \item{\code{call} : }{The output of \code{\link{match.call}} for the user-
#'      supplied function arguments.}
#' }
#'
#' @details This function is a \code{\link{switch}} call to
#'    \code{\link{coxTrain_fun}} (for \code{type = "survival"}),
#'    \code{\link{olsTrain_fun}} (for \code{type = "regression"}), or
#'    \code{\link{glmTrain_fun}} (for \code{type = "categorical"}).
#'
#' @seealso \code{\link{superpc.st}}; \code{\link{superPCA_pVals}}
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use superPCA_pVals() instead

superpc.train <- function(data,
                          type = c("survival", "regression", "categorical"),
                          s0.perc = NULL){

  # browser()

  this.call <- match.call()
  type <- match.arg(type)

  ###  Error Checks  ###
  censor_logi <- is.null(data$censoring.status)

  if(censor_logi & type == "survival"){
    stop("Error: survival specified but censoring.status is null")
  }

  if(!censor_logi & type == "regression"){
    stop("Error: regression specified but censoring.status is non-null")
  }

  ###  Model Switch  ###
  switch(type,
    survival = {
      junk <- coxTrain_fun(
        data$x, data$y,
        data$censoring.status,
        s0.perc = s0.perc
      )
    },
    regression = {
      junk <- olsTrain_fun(
        as.matrix(data$x), data$y,
        s0.perc = s0.perc
      )
    },
    categorical = {

      resp <- data$y

      if(!(is.integer(resp) | is.factor(resp))){
        stop("Response must be an integer or factor for classification.")
      }

      if(is.ordered(resp)){

        stop("Ordered Logistic Regression not currently implemented.")
        type <- "ordered"
        # MASS::polr implementation

      } else if(length(unique(resp)) > 2) {

        stop("Multinomial Regression not currently implemented.")
        type <- "n_ary"
        # nnet::multinom implementation

      } else if(length(unique(resp)) == 2) {

        type <- "binary"
        junk <- glmTrain_fun(data$x, resp, family = binomial)
        junk$tt[is.na(junk$tt)] <- 0

      } else {
        stop("Response requires two or more distinct values for classification.")
      }

    })

  out_ls <- list(feature.scores = junk$tt,
                 type = type,
                 s0.perc = s0.perc,
                 call = this.call)


  class(out_ls) <- "superpc"
  out_ls

}

