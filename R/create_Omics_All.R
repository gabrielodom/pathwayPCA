#' Generation Functions for -Omics.* Class Objects
#'
#' These functions create valid objects of the "OmicsPathway", "OmicsSurv",
#'   "OmicsReg", and "OmicsCateg" classes.
#'
#' @details Please note that the classes of the parameters are \emph{not}
#'   flexible. The -Omics measurement data \emph{must} be a data frame, and the
#'   response values (for the survival, regression, and classification object)
#'   \emph{must} match their expected classes \emph{exactly}. The reason for
#'   this is to force the end user to pay attention to the quality and format
#'   of their input data. Because the functions internal to this package have
#'   only been tested on the classes described in the Arguments section, these
#'   class checks prevent unexpected errors (or worse, incorrect computational
#'   results without an error). These draconian input class restrictions
#'   protect the accuracy of your data's analysis.
#'
#' @section OmicsPathway:
#' Valid OmicsPathawy objects will have no response information, just the mass
#'   spectrometry (gene "design") matrix and the pathway list. OmicsPathway
#'   objects should be created only when unsupervised pathway extraction is
#'   needed. Because of the missing response, no pathway testing can be
#'   performed on an OmicsPathway object.
#'
#' @param assayData_df An $N x p$ data frame with named columns
#' @param pathwaySet_ls A list of known gene pathways with one or two elements:
#' \itemize{
#'   \item{pathways : }{A named list of character vectors. Each vector contains
#'     the names of the individual genes within that pathway as a vector of
#'     character strings. The names contained in these vectors must have non-
#'     empty overlap with the \emph{column names} of the \code{assayData_df}
#'     data frame. The names of the pathways (the list elements themselves)
#'     should be the a shorthand representation of the full pathway name.}
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
#'   colon_OmicsPath <- create_OmicsPath(assayData_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colonGenesets_ls)
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_OmicsSurv(assayData_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colonGenesets_ls,
#'                                       eventTime_vec = colonSurv_df$OS_time,
#'                                       eventObserved_vec = as.logical(colonSurv_df$OS_event))
#'
#'   ###  Create an OmicsReg Object  ###
#'   colon_OmicsReg <- create_OmicsReg(assayData_df = colonSurv_df[, -(1:2)],
#'                                     pathwaySet_ls = colonGenesets_ls,
#'                                     response_num = colonSurv_df$OS_time)
#'
#'   ###  Create an OmicsCateg Object  ###
#'   colon_OmicsCateg <- create_OmicsCateg(assayData_df = colonSurv_df[, -(1:2)],
#'                                         pathwaySet_ls = colonGenesets_ls,
#'                                         response_fact = as.factor(colonSurv_df$OS_event))
#'
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsPath <- function(assayData_df, pathwaySet_ls){

  if("matrix" %in% class(assayData_df) &
     !("data.frame" %in% class(assayData_df))){
    stop("\n You have supplied a matrix object to the assayData_df argument.
    Note that the pathwayPCA:: package functions require -Omics data as an N x p
    data frame object: this data frame will have one observation per row and one
    measurement per column. If your matrix is in 'tall' (p x N) format, please
    transpose your matrix with the 't()' function (but pay attention to your
    column names after transposition). Next, you can use the 'as.data.frame()'
    function to transform your -Omics data matrix to class 'data.frame'. Please
    see the help information found in ?create_OmicsPath for more details.")
  }

  pathwaySet_ls$setsize <- sapply(pathwaySet_ls$pathways, length)

  new("OmicsPathway",
      assayData_df = assayData_df,
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
create_OmicsSurv <- function(assayData_df,
                             pathwaySet_ls,
                             eventTime_vec,
                             eventObserved_vec){

  if("matrix" %in% class(assayData_df) &
     !("data.frame" %in% class(assayData_df))){
    stop("\n You have supplied a matrix object to the assayData_df argument.
    Note that the pathwayPCA:: package functions require -Omics data as an N x p
    data frame object: this data frame will have one observation per row and one
    measurement per column. If your matrix is in 'tall' (p x N) format, please
    transpose your matrix with the 't()' function (but pay attention to your
    column names after transposition). Next, you can use the 'as.data.frame()'
    function to transform your -Omics data matrix to class 'data.frame'. Please
    see the help information found in ?create_OmicsSurv for more details.")
  }

  pathwaySet_ls$setsize <- sapply(pathwaySet_ls$pathways, length)

  new("OmicsSurv",
      assayData_df = assayData_df,
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
create_OmicsReg <- function(assayData_df,
                            pathwaySet_ls,
                            response_num){

  if("matrix" %in% class(assayData_df) &
     !("data.frame" %in% class(assayData_df))){
    stop("\n You have supplied a matrix object to the assayData_df argument.
    Note that the pathwayPCA:: package functions require -Omics data as an N x p
    data frame object: this data frame will have one observation per row and one
    measurement per column. If your matrix is in 'tall' (p x N) format, please
    transpose your matrix with the 't()' function (but pay attention to your
    column names after transposition). Next, you can use the 'as.data.frame()'
    function to transform your -Omics data matrix to class 'data.frame'. Please
    see the help information found in ?create_OmicsReg for more details.")
  }

  pathwaySet_ls$setsize <- sapply(pathwaySet_ls$pathways, length)

  new("OmicsReg",
      assayData_df = assayData_df,
      pathwaySet = pathwaySet_ls,
      response = response_num)

}

#' @param response_fact A factor vector of length $N$: the dependent variable of a
#'   generalized linear model
#' @export
#' @rdname create_OmicsPathway
create_OmicsCateg <- function(assayData_df,
                              pathwaySet_ls,
                              response_fact){

  if("matrix" %in% class(assayData_df) &
     !("data.frame" %in% class(assayData_df))){
    stop("\n You have supplied a matrix object to the assayData_df argument.
    Note that the pathwayPCA:: package functions require -Omics data as an N x p
    data frame object: this data frame will have one observation per row and one
    measurement per column. If your matrix is in 'tall' (p x N) format, please
    transpose your matrix with the 't()' function (but pay attention to your
    column names after transposition). Next, you can use the 'as.data.frame()'
    function to transform your -Omics data matrix to class 'data.frame'. Please
    see the help information found in ?create_OmicsCateg for more details.")
  }

  pathwaySet_ls$setsize <- sapply(pathwaySet_ls$pathways, length)

  new("OmicsCateg",
      assayData_df = assayData_df,
      pathwaySet = pathwaySet_ls,
      response = response_fact)

}
