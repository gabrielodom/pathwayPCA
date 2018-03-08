#' Generation Functions for -Omics.* Class Objects
#'
#' These functions create valid objects of the "OmicsPathway", "OmicsSurv",
#'   "OmicsReg", and "OmicsCateg" classes.
#'
#' @section OmicsPathway:
#' Valid OmicsPathawy objects will have no response information, just the mass
#'   spectrometry (gene "design") matrix and the pathway list. OmicsPathway
#'   objects should be created only when unsupervised pathway extraction is
#'   needed. Because of the missing response, no pathway testing can be
#'   performed on an OmicsPathway object.
#'
#' @param massSpec_df An $N x p$ data frame with named columns
#' @param pathwaySet_ls A list of known gene pathways with one or two elements:
#' \itemize{
#'   \item{pathways : }{A named list of character vectors. Each vector contains
#'     the names of the individual genes within that pathway as a vector of
#'     character strings. The names contained in these vectors must have non-
#'     empty overlap with the \emph{column names} of the \code{massSpec} data
#'     frame. The names of the pathways (the list elements themselves) should
#'     be the a shorthand representation of the full pathway name.}
#'   \item{TERMS: }{ A character vector the same length as the
#'     \code{pathways} list with the proper names of the pathways.}
#' }
#'
#' @return An object of classes "OmicsPathway", "OmicsSurv", "OmicsReg", or
#'   "OmicsCateg".
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#' @include createClass_OmicsSurv.R
#' @include createClass_OmicsReg.R
#' @include createClass_OmicsCateg.R
#'
#' @seealso \code{"\link[=OmicsPathway-class]{OmicsPathway}"},
#'   \code{"\link[=OmicsSurv-class]{OmicsSurv}"},
#'   \code{"\link[=OmicsReg-class]{OmicsReg}"},  and
#'   \code{"\link[=OmicsCateg-class]{OmicsCateg}"}
#'
#' @importFrom methods new
#'
#' @examples
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colonGenesets_ls")
#'
#'   ###  Create an OmicsPathway Object  ###
#'   colon_OmicsPath <- create_OmicsPath(massSpec_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colonGenesets_ls)
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_OmicsSurv(massSpec_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colonGenesets_ls,
#'                                       eventTime_vec = colonSurv_df$OS_time,
#'                                       eventObserved_vec = as.logical(colonSurv_df$OS_event))
#'
#'   ###  Create an OmicsReg Object  ###
#'   colon_OmicsReg <- create_OmicsReg(massSpec_df = colonSurv_df[, -(1:2)],
#'                                     pathwaySet_ls = colonGenesets_ls,
#'                                     response_num = colonSurv_df$OS_time)
#'
#'   ###  Create an OmicsCateg Object  ###
#'   colon_OmicsCateg <- create_OmicsCateg(massSpec_df = colonSurv_df[, -(1:2)],
#'                                         pathwaySet_ls = colonGenesets_ls,
#'                                         response_fact = as.factor(colonSurv_df$OS_event))
#'
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsPath <- function(massSpec_df, pathwaySet_ls){

  pathwaySet_ls$setsize <- unname(sapply(pathwaySet_ls$pathways, length))

  new("OmicsPathway",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls)

}

#' @section OmicsSurv:
#' Valid OmicsSurv objects will have two response vectors: a vector of the most
#'   recently recorded follow-up times and a logical vector if that time marks
#'   an event (TRUE = observed event; FALSE = right-censored observation).
#'
#' @param eventTime_vec A numeric vector with $N$ observations corresponding to
#'   the last observed time of follow up
#' @param eventObserved_vec A logical vector with $N$ observations indicating
#'   right censoring. The values will be FALSE if the observation was censored
#'   (i.e., we did not observe an event).
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsSurv <- function(massSpec_df,
                             pathwaySet_ls,
                             eventTime_vec,
                             eventObserved_vec){

  pathwaySet_ls$setsize <- unname(sapply(pathwaySet_ls$pathways, length))

  new("OmicsSurv",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls,
      eventTime = eventTime_vec,
      eventObserved = eventObserved_vec)

}




#' @section OmicsReg and OmicsCateg:
#' Valid OmicsReg and OmicsCateg objects with have one response vector of
#'   continuous or categorial (as a factor) observations, respectively.
#'
#' @param response_num A numeric vector of length $N$: the dependent variable
#'   in a regression exercise
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsReg <- function(massSpec_df,
                            pathwaySet_ls,
                            response_num){

  pathwaySet_ls$setsize <- unname(sapply(pathwaySet_ls$pathways, length))

  new("OmicsReg",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls,
      response = response_num)

}

#' @param response_fact A factor vector of length $N$: the dependent variable of a
#'   generalized linear model
#' @export
#' @rdname create_OmicsPathway
create_OmicsCateg <- function(massSpec_df,
                              pathwaySet_ls,
                              response_fact){

  pathwaySet_ls$setsize <- unname(sapply(pathwaySet_ls$pathways, length))

  new("OmicsCateg",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls,
      response = response_fact)

}
