#' Validity Checking for -Omics.* Classes
#'
#' These functions check the validity of the "OmicsSurv", "OmicsReg", and
#'   "OmicsCateg" classes.
#'
#' At the moment, we have currently written checks to make sure the dimensions
#'   of the mass spectrometry or bio-assay data frame and response vectors
#'   match.
#'
#' @section OmicsSurv:
#' Valid OmicsSurv objects will have two response vectors: a vector of the most
#'   recently recorded follow-up times and a logical vector if that time marks
#'   an event (TRUE = observed event; FALSE = right-censored observation).
#'
#' @param object An object of classes "OmicsSurv", "OmicsReg", or "OmicsCateg".
#' @return TRUE if the object is a valid object, else an error message with the
#'   rule broken.
#' @rdname valid_OmicsSurv
valid_OmicsSurv <- function(object){

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
#' Valid OmicsReg and OmicsCateg objects with have one response vector of
#'   continuous or categorial (as a factor) observations, respectively.
#' @rdname valid_OmicsSurv
valid_OmicsReg <- function(object){

  nX <- nrow(object@assayData_df)
  nY <- length(object@response)

  if(nY != nX){
    return("Number of assayData_df rows must match number of responses.")
  } else {
    return(TRUE)
  }

}

#' @rdname valid_OmicsSurv
valid_OmicsCateg <- function(object){

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



# library(methods)
# library(pathwayPCA)
# load("data/ovarianFiltered_df.rda")
# load("data/genesets_ls.rda")
