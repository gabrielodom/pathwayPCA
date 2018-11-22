#' Extract PCs and Loadings from a \code{superpcOut}- or \code{aespcOut}-class
#'    Object.
#'
#' @description Given an object of class \code{aespcOut} or \code{superpcOut},
#'    as returned by the functions \code{\link{AESPCA_pVals}} or
#'    \code{\link{SuperPCA_pVals}}, respectively, and the name or unique ID of
#'    a pathway, return a data frame of the principal components and a data
#'    frame of the loading vectors corresponding to that pathway.
#'
#' @details Match the supplied pathway character string to either the
#'    \code{pathways} or \code{terms} columns of the \code{pVals_df} data frame
#'    within the \code{pcOut} object. Then, subset the \code{loadings_ls} and
#'    \code{PCs_ls} lists for their entries which match the supplied pathway.
#'    Finally, return a list of the PCs, loadings, and the pathway ID and name.
#'
#' @param pcOut An object of classes \code{superpcOut} or \code{aespcOut} as
#'    returned by the \code{\link{SuperPCA_pVals}} or \code{\link{AESPCA_pVals}}
#'    functions, respectively.
#' @param pathway_char A character string of the name or unique identifier of a
#'    pathway
#' @param ... Dots for additional arguments (currently unused).
#'
#' @return A list of four elements:
#'    \itemize{
#'      \item{\code{PCs} : }{A data frame of the principal components}
#'      \item{\code{Loadings} : }{A matrix of the loading vectors with features
#'        in the row names}
#'      \item{\code{pathway} : }{The unique pathway identifier for the
#'        \code{pcOut} object}
#'      \item{\code{term} : }{The name of the pathway}
#'    }
#'
#' @examples
#'   NULL
#'
#' @rdname getPathSVD
#' @export getPathSVD
getPathSVD <- function(pcOut, pathway_char, ...){
  UseMethod("getPathSVD")
}

#' @return \code{NULL}
#'
#' @rdname getPathSVD
#' @method getPathSVD superpcOut
#' @S3method getPathSVD superpcOut
getPathSVD.superpcOut <- function(pcOut, pathway_char, ...){

  ###  Check for Matches  ###
  pathID_idx <- which(pcOut$pVals_df$pathways == pathway_char)
  term_idx   <- which(pcOut$pVals_df$terms == pathway_char)
  if(length(pathID_idx) + length(term_idx) == 0){
    stop("Supplied pathway does not match any pathway in the supplied object.")
  }

  ###  Match pathway ID to pathway name  ###
  if(length(term_idx) == 1){

    pathID  <- as.character(pcOut$pVals_df[term_idx, "pathways"])
    pathway <- pathway_char

  } else if(length(pathID_idx) == 1) {

    pathID  <- pathway_char
    pathway <- as.character(pcOut$pVals_df[pathID_idx, "terms"])

  }

  ###  Select the PCs and Loadings that Match this Pathway ID  ###
  list(
    PCs = pcOut$PCs_ls[[pathID]],
    Loadings = pcOut$loadings_ls[[pathID]],
    pathway = pathID,
    term = pathway
  )

}

#' @return \code{NULL}
#'
#' @rdname getPathSVD
#' @method getPathSVD aespcOut
#' @S3method getPathSVD aespcOut
getPathSVD.aespcOut <- getPathSVD.superpcOut
