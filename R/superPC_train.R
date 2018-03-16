#' Train Supervised PC Model
#'
#' @description Computes feature scores for supervised pc analysis
#'
#' @param data A list of training data:
#' \itemize{
#'   \item{x : }{A "tall" pathway data frame ($p_{path} * n$).}
#'   \item{y : }{A response vector corresponding to \code{type}}
#'   \item{censoring.status : }{If \code{type = "survival"}, the censoring
#'     indicator. Otherwise, \code{NULL}.}
#'   \item{featurenames : }{A character vector of the measured -omes in
#'     \code{x}.}
#' }
#' @param type What model relates \code{y} and \code{x}? Options are
#'    \code{"survival"}, \code{"regression"}, or \code{"classification"} (for
#'    (potentially multinomial)logistic  regression).
#' @param s0.perc A stabilization parameter on the interval [0,1]. This is an
#'    internal argument to each of the called functions. The default NULL value
#'    will ensure an appropriate value is determined internally.
#'
#' @return A list containing:
#' \itemize{
#'   \item{feature.scores : }{The scaled p-dimensional score vector: each value
#'     has been divided by its respective standard deviation plus epsilon
#'     (governed by \code{s0.perc}). \code{NA} values returned by the logistic
#'     model are replaced with 0.}
#'   \item{type : }{The argument for \code{type}.}
#'   \item{s0.perc : }{The user-supplied value of \code{s0.perc}, or the
#'     internally-calculated default value from the chosen model.}
#'   \item{call : }{The output of \code{match.call()}.}
#' }
#'
#' @details This function is a switch call to \code{\link{coxTrain_fun}},
#'    \code{\link{olsTrain_fun}}, or \code{glmTrain_fun}, respectively.
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use superPCA_pVals() instead

superpc.train <- function(data,
                          type = c("survival", "regression", "classification"),
                          s0.perc = NULL){

  # browser()

  this.call <- match.call()
  type <- match.arg(type)

  ###  Error Checks  ###
  censor_logic <- is.null(data$censoring.status)

  if(censor_logic & type == "survival"){
    stop("Error: survival specified but censoring.status is null")
  }

  if(!censor_logic & type == "regression"){
    stop("Error: regression specified but censoring.status is non-null")
  }

  ###  Model Switch  ###
  switch(type,
         survival = {

           junk <- coxTrain_fun(data$x, data$y, data$censoring.status,
                                s0.perc = s0.perc)
           feature.scores <- junk$tt

         },
         regression = {

           junk <- olsTrain_fun(as.matrix(data$x), data$y, s0.perc = s0.perc)
           feature.scores <- junk$tt

         },
         classification = {

           resp <- data$y
           if(!(is.integer(resp) | is.factor(resp))){
             stop("Response must be an integer or factor for classification.")
           }

           if(is.ordered(resp)){

             stop("Ordered Logistic Regression not currently implemented.")
             type <- "ordered"
             # MASS::polr implementation

           } else if(length(unique(resp)) > 2){

             stop("Multinomial Regression not currently implemented.")
             type <- "n_ary"
             # nnet::multinom implementation

           } else if(length(unique(resp)) == 2){

             type <- "binary"
             junk <- glmTrain_fun(data$x, resp, family = binomial)
             feature.scores <- junk$tt
             feature.scores[is.na(feature.scores)] <- 0

           } else {
             stop("Response requires two or more distinct values for classification.")
           }



         })

  junk <- list(feature.scores = feature.scores,
               type = type,
               s0.perc = s0.perc,
               call = this.call)


  class(junk) <- "superpc"
  return(junk)

}

