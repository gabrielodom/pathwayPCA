#' Extract Table of \eqn{p}-values from a \code{superpcOut}- or \code{aespcOut}-
#'    class Object.
#'
#' @description Given an object of class \code{aespcOut} or \code{superpcOut},
#'    as returned by the functions \code{\link{AESPCA_pVals}} or
#'    \code{\link{SuperPCA_pVals}}, respectively, return a data frame of the
#'    \eqn{p}-values for the top pathways.
#'
#' @details Row-subset the \code{pVals_df} entry of an object of class
#'    \code{aespcOut} or \code{superpcOut} by the number of pathways requested
#'    (via the \code{nPaths} argument) or by the unadjusted significance level
#'    for each pathway (via the \code{alpha} argument). Return a data frame of
#'    the pathway names, FDR-adjusted significance levels (if available), and
#'    the raw score (negative natural logarithm of the \eqn{p}-values) of each
#'    pathway.
#'
#' @param pcOut An object of classes \code{superpcOut} or \code{aespcOut} as
#'    returned by the \code{\link{SuperPCA_pVals}} or \code{\link{AESPCA_pVals}}
#'    functions, respectively.
#' @param score Should the unadjusted \eqn{p}-values be returned transformed to
#'    negative natural logarithm scores or left as is? Defaults to \code{FALSE};
#'    that is, the raw \eqn{p}-values are returned instead of the transformed
#'    \eqn{p}-values.
#' @param numPaths The number of top pathways by raw \eqn{p}-value. Defaults to
#'    the top 20 pathways. We do not permit users to specify \code{numPaths}
#'    and \code{alpha} concurrently.
#' @param alpha The significance threshold for raw \eqn{p}-values. Defaults to
#'    \code{NULL}. If \code{alpha} is given, then \code{numPaths} will be
#'    ignored.
#' @param ... Dots for additional arguments (currently unused).
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item{\code{terms} : }{The pathway name, as given in the
#'       \code{object@@trimPathwayCollection$TERMS} object.}
#'     \item{\code{description} : }{(OPTIONAL) The pathway description, as given
#'       in the \code{object@@trimPathwayCollection$description} object, if
#'       supplied.}
#'     \item{\code{rawp} : }{The unadjusted \eqn{p}-values of each pathway.
#'       Included if \code{score = FALSE}.}
#'     \item{\code{...} : }{Additional columns of FDR-adjusted \eqn{p}-values
#'       as specified through the \code{adjustment} argument of the
#'       \code{\link{SuperPCA_pVals}} or \code{\link{AESPCA_pVals}} functions.}
#'     \item{\code{score} : }{The negative natural logarithm of the unadjusted
#'       \eqn{p}-values of each pathway. Included if \code{score = TRUE}.}
#'   }
#'
#' @examples
#'
#'   ###  Load Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create -Omics Container  ###
#'   colon_Omics <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "survival"
#'   )
#'
#'   ###  Calculate Supervised PCA Pathway p-Values  ###
#'   colon_superpc <- SuperPCA_pVals(
#'     colon_Omics,
#'     numPCs = 2,
#'     parallel = TRUE,
#'     numCores = 2,
#'     adjustment = "BH"
#'   )
#'
#'   ###  Extract Table of p-Values  ###
#'   # Top 5 Pathways
#'   getPathpVals(
#'     colon_superpc,
#'     numPaths = 5
#'   )
#'   
#'   # Pathways with Unadjusted p-Values < 0.01
#'   getPathpVals(
#'     colon_superpc,
#'     alpha = 0.01
#'   )
#'
#'
#' @rdname getPathpVals
#' @export getPathpVals
getPathpVals <- function(pcOut, score = FALSE,
                         numPaths = 20L, alpha = NULL, ...){
  UseMethod("getPathpVals")
}

#' @return \code{NULL}
#'
#' @rdname getPathpVals
#' @method getPathpVals superpcOut
#' @export
getPathpVals.superpcOut <- function(pcOut, score = FALSE,
                                    numPaths = 20L, alpha = NULL, ...){
  
  pVals_df <- pcOut$pVals_df
  # The output from AESPCA_pVals() and SuperPCA_pVals() always contains a 
  #   *sorted* data frame of p-values, so we don't have to do any sorting here.
  
  ###  Select the Requested Rows  ###
  if(is.null(alpha)){
    
    numPaths <- min(nrow(pVals_df), numPaths)
    pValsShort_df <- pVals_df[seq_len(numPaths), ]
    
  } else {
    
    if(!missing(numPaths)){
      warning(
        "numPaths and alpha cannot both be specified. Using alpha.",
        call. = FALSE
      )
    }
    
    pValsShort_df <- pVals_df[pVals_df$rawp < alpha, ]
    
    if(nrow(pValsShort_df) == 0){
      return(pValsShort_df)
    }
    
  }
  
  ###  Add / Delete Variables  ###
  if(score){
    
    pValsShort_df$score  <- -log(pValsShort_df$rawp)
    pValsShort_df$rawp   <- NULL
    
  }
  
  pValsShort_df$pathways <- NULL
  pValsShort_df$n_tested <- NULL
  
  
  ###  Return  ###
  pValsShort_df
  
}

#' @return \code{NULL}
#'
#' @rdname getPathpVals
#' @method getPathpVals aespcOut
#' @export
getPathpVals.aespcOut <- getPathpVals.superpcOut
