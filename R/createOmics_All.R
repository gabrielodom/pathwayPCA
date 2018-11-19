#' Generation functions for \code{-Omics*}-class objects
#'
#' These functions create valid objects of class \code{OmicsPathway},
#'   \code{OmicsSurv}, \code{OmicsReg}, or \code{OmicsCateg}.
#'
#' @details Please note that the classes of the parameters are \emph{not}
#'   flexible. The -Omics assay data \emph{must} be or extend the class
#'   \code{data.frame}, and the response values (for a survival-, regression-,
#'   or categorical-response object) \emph{must} match their expected classes
#'   \emph{exactly}. The reason for this is to encourage the end user to pay
#'   attention to the quality and format of their input data. Because the
#'   functions internal to this package have only been tested on the classes
#'   described in the Arguments section, these class checks prevent unexpected
#'   errors (or worse, incorrect computational results without an error). These
#'   draconian input class restrictions protect the accuracy of your data
#'   analysis.
#'
#' @section OmicsPathway:
#' Valid \code{OmicsPathway} objects will have no response information, just the
#'   mass spectrometry or bio-assay ("design") matrix and the pathway list.
#'   \code{OmicsPathway} objects should be created only when unsupervised
#'   pathway extraction is needed (not possible with Supervised PCA). Because of
#'   the missing response, no pathway testing can be performed on an
#'   \code{OmicsPathway} object.
#'
#' @param assayData_df An \eqn{N \times p} data frame with named columns.
#' @param pathwayCollection_ls A \code{pathwayCollection} list of known gene
#'   pathways with two or three elements:
#'   \itemize{
#'     \item{\code{pathways} : }{A named list of character vectors. Each vector
#'        contains the names of the individual genes within that pathway as a
#'        vector of character strings. The names contained in these vectors must
#'        have non-empty overlap with the \emph{column names} of the
#'        \code{assayData_df} data frame. The names of the pathways (the list
#'        elements themselves) should be the a shorthand representation of the
#'        full pathway name.}
#'     \item{\code{TERMS}: }{A character vector the same length as the
#'        \code{pathways} list with the proper names of the pathways.}
#'     \item{\code{description} : }{An optional character vector the same length
#'        as the \code{pathways} list with additional information about the
#'        pathways.}
#'   }
#'
#' @return A valid object of class \code{OmicsPathway}, \code{OmicsSurv},
#'   \code{OmicsReg}, or \code{OmicsCateg}.
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#' @include createClass_OmicsSurv.R
#' @include createClass_OmicsReg.R
#' @include createClass_OmicsCateg.R
#'
#' @seealso \code{\link[=OmicsPathway-class]{OmicsPathway}},
#'   \code{\link[=OmicsSurv-class]{OmicsSurv}},
#'   \code{\link[=OmicsReg-class]{OmicsReg}}, and
#'   \code{\link[=OmicsCateg-class]{OmicsCateg}}
#'
#' @importFrom methods new
#'
#' @examples
#'   DO NOT CALL THESE FUNCTIONS DIRECTLY. USE CreateOmics() INSTEAD.
#'
#' @export
#' @rdname CreateOmicsPathway
CreateOmicsPath <- function(assayData_df,
                            pathwayCollection_ls){

  obj <- new(
    "OmicsPathway",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls
  )

}




#' @section OmicsSurv:
#' Valid \code{OmicsSurv} objects will have two response vectors: a vector of
#'   the most recently recorded follow-up times and a logical vector if that
#'   time marks an event (\code{TRUE}: observed event; \code{FALSE}: right-
#'   censored observation).
#'
#' @param eventTime_num A \code{numeric} vector with \eqn{N} observations
#'    corresponding to the last observed time of follow up.
#' @param eventObserved_lgl A \code{logical} vector with \eqn{N} observations
#'   indicating right-censoring. The values will be \code{FALSE} if the
#'   observation was censored (i.e., we did not observe an event).
#'
#' @export
#' @rdname CreateOmicsPathway
CreateOmicsSurv <- function(assayData_df,
                            pathwayCollection_ls,
                            eventTime_num,
                            eventObserved_lgl){

  obj <- new(
    "OmicsSurv",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls,
    eventTime = eventTime_num,
    eventObserved = eventObserved_lgl
  )

}




#' @section OmicsReg and OmicsCateg:
#' Valid \code{OmicsReg} and \code{OmicsCateg} objects with have one response
#'   vector of continuous (\code{numeric}) or categorial (\code{factor})
#'   observations, respectively.
#'
#' @param response_num A \code{numeric} vector of length \eqn{N}: the dependent
#'   variable in an ordinary regression exercise.
#'
#' @export
#' @rdname CreateOmicsPathway
CreateOmicsReg <- function(assayData_df,
                            pathwayCollection_ls,
                            response_num){

  obj <- new(
    "OmicsReg",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls,
    response = response_num
  )

}




#' @param response_fact A \code{factor} vector of length \eqn{N}: the dependent
#'   variable of a generalized linear regression exercise.
#' @export
#' @rdname CreateOmicsPathway
CreateOmicsCateg <- function(assayData_df,
                              pathwayCollection_ls,
                              response_fact){

  obj <- new(
    "OmicsCateg",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls,
    response = response_fact
  )

}
