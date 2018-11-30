#' Check validity of new Omics*-class objects
#'
#' @description These functions check the validity of new objects created in
#'    the \code{OmicsSurv}, \code{OmicsReg}, and \code{OmicsCateg} classes.
#'
#' @details We have currently written checks to make sure the dimensions of the
#'    mass spectrometry or bio-assay data frame and response matrices or vectors
#'    match. Other checks should be added in response to user feedback during or
#'    after beta testing. ENHANCEMENT.
#'
#' @keywords internal
#'
#' @section OmicsSurv:
#' Valid \code{OmicsSurv} objects will have two response vectors: a vector of
#'    the most recently recorded follow-up times and a logical vector if that
#'    time marks a death or event (\code{TRUE}: observed event; \code{FALSE}:
#'    right-censored observation).
#'
#' @param object An object potentially of class \code{OmicsSurv},
#'    \code{OmicsReg}, or \code{OmicsCateg}.
#'
#' @return \code{TRUE} if the object is a valid object, else an error message
#'    with the rule broken.
#'
#' @rdname ValidOmicsSurv
ValidOmicsSurv <- function(object){

  nX <- nrow(object@assayData_df)
  nY <- length(object@eventTime)
  nDelta <- length(object@eventObserved)

  if(nY != nX){
    return("Number of assayData_df rows must match number of response times.")
  } else if(nDelta != nX){
    return("Number of assayData_df rows must match number of response censoring indicators.")
  } else {
    return(TRUE)
  }

}

#' @section OmicsReg and OmicsCateg:
#' Valid \code{OmicsReg} and \code{OmicsCateg} objects with have one response
#'    vector of continuous (\code{numeric}) or categorial (\code{factor})
#'    observations, respectively.
#'
#' @rdname ValidOmicsSurv
ValidOmicsReg <- function(object){

  nX <- nrow(object@assayData_df)
  nY <- length(object@response)

  if(nY != nX){
    return("Number of assayData_df rows must match number of responses.")
  } else {
    return(TRUE)
  }

}

#' @rdname ValidOmicsSurv
ValidOmicsCateg <- function(object){

  nX <- nrow(object@assayData_df)
  nY <- length(object@response)

  if(nlevels(object@response) > 2){
    return("\n The current implementation of the pathwayPCA package permits binary response only.")
  }

  if(nY != nX){
    return("Number of assayData_df rows must match number of responses.")
  } else {
    return(TRUE)
  }

}
