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
#'   Some of the pathways in the supplied pathways list will be removed, or
#'   "trimmed", during object creation. For the pathway-testing methods, these
#'   trimmed pathways will have \eqn{p}-values given as \code{NA}. For an
#'   explanation of pathway trimming, see the documentation for the
#'   \code{\link{IntersectOmicsPwyCollct}} function.
#'
#'   Also note the following: if the supplied \code{pathways} object within your
#'   \code{pathwayCollection_ls} list has no names, then this pathway list will
#'   be named \code{path1}, \code{path2}, \code{path3}, ...; if any of the
#'   pathways are missing names, then the missing pathways will be named
#'   \code{noName} followed by the index of the pathway. For example, if the
#'   112th pathway in the \code{pathways} list has no name (but other pathways
#'   do), then this pathway will be named \code{noName112}. Furthermore, if any
#'   of the pathway names are duplicated, then the duplicates will have
#'   \code{.1}, \code{.2}, \code{.3}, ... appended to the duplicate names until
#'   all pathway names are unique. Once all pathways have been verified to have
#'   unique names, then the pathway names are attached as attributes to the
#'   \code{TERMS} and \code{setsize} vectors (the \code{setsize} vector is
#'   calculated at object creation).
#'
#'   If your gene pathways list is stored in a \code{.gmt} file, use the
#'   \code{\link{read_gmt}} function to import your pathways list as a
#'   \code{pathwayCollection} list object.
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
#' @param minPathSize What is the smallest number of genes allowed in each
#'   pathway? Defaults to 3.
#' @param ... Dots for additional arguments passed to the internal
#'   \code{\link{CheckAssay}} function.
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
#'   \code{\link[=OmicsReg-class]{OmicsReg}},
#'   \code{\link[=OmicsCateg-class]{OmicsCateg}}, and
#'   \code{\link{IntersectOmicsPwyCollct}}
#'
#' @importFrom methods new
#'
#' @examples
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create an OmicsPathway Object  ###
#'   colon_OmicsPath <- CreateOmicsPath(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection
#'   )
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- CreateOmicsSurv(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     eventTime_num = colonSurv_df$OS_time,
#'     eventObserved_lgl = as.logical(colonSurv_df$OS_event)
#'   )
#'
#'   ###  Create an OmicsReg Object  ###
#'   colon_OmicsReg <- CreateOmicsReg(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response_num = colonSurv_df$OS_time
#'   )
#'
#'   ###  Create an OmicsCateg Object  ###
#'   colon_OmicsCateg <- CreateOmicsCateg(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response_fact = as.factor(colonSurv_df$OS_event)
#'   )
#' }
#'
#' @export
#' @rdname CreateOmicsPathway
CreateOmicsPath <- function(assayData_df,
                            pathwayCollection_ls,
                            minPathSize = 3,
                            ...){
  # browser()

  ###  Error Checks and Warnings for supplied assay  ###
  assayData_df <- CheckAssay(assayData_df, ...)


  ###  Process Pathway Collection  ###
  pathwayCollection_ls$setsize <- lengths(pathwayCollection_ls$pathways)

  # If there are no names, create them. If there are missing names, label them.
  if(is.null(names(pathwayCollection_ls$pathways))){

    pathNames <- paste0("path", 1:length(pathwayCollection_ls$pathways))
    names(pathwayCollection_ls$pathways) <- pathNames

  } else if(anyNA(names(pathwayCollection_ls$pathways))){

    missingNm_idx <- which(is.na(names(pathwayCollection_ls$pathways)))
    names(pathwayCollection_ls$pathways)[missingNm_idx] <-
      paste0("noName", missingNm_idx)
    pathNames <- names(pathwayCollection_ls$pathways)

  } else {
    pathNames <- names(pathwayCollection_ls$pathways)
  }

  if(anyDuplicated(names(pathwayCollection_ls$pathways))){

    shortSet_df <- data.frame(
      lapply(pathwayCollection_ls$pathways, "length<-", 1)
    )
    pathNames <- names(shortSet_df)
    names(pathwayCollection_ls$pathways) <- pathNames

  }

  # Add Name Key to TERMS and setsize
  names(pathwayCollection_ls$TERMS) <- pathNames
  names(pathwayCollection_ls$setsize) <- pathNames


  ###  Create Omics Object  ###
  obj <- new(
    "OmicsPathway",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls
  )

  IntersectOmicsPwyCollct(obj, trim = minPathSize)

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
                            eventObserved_lgl,
                            minPathSize = 3,
                            ...){

  ###  Error Checks and Warnings for supplied assay  ###
  assayData_df <- CheckAssay(assayData_df, ...)


  ###  Process Pathway Collection  ###
  pathwayCollection_ls$setsize <- lengths(pathwayCollection_ls$pathways)


  ###  Pathways List Names Checking  ###
  # If there are no names, create them. If there are missing names, label them.
  if(is.null(names(pathwayCollection_ls$pathways))){

    pathNames <- paste0("path", 1:length(pathwayCollection_ls$pathways))
    names(pathwayCollection_ls$pathways) <- pathNames

  } else if(anyNA(names(pathwayCollection_ls$pathways))){

    missingNm_idx <- which(is.na(names(pathwayCollection_ls$pathways)))
    names(pathwayCollection_ls$pathways)[missingNm_idx] <-
      paste0("noName", missingNm_idx)
    pathNames <- names(pathwayCollection_ls$pathways)

  } else {
    pathNames <- names(pathwayCollection_ls$pathways)
  }

  if(anyDuplicated(names(pathwayCollection_ls$pathways))){

    shortSet_df <- data.frame(
      lapply(pathwayCollection_ls$pathways, "length<-", 1)
    )
    pathNames <- names(shortSet_df)
    names(pathwayCollection_ls$pathways) <- pathNames

  }


  ###  Add Name Key to TERMS and setsize  ###
  names(pathwayCollection_ls$TERMS) <- pathNames
  names(pathwayCollection_ls$setsize) <- pathNames


  ###  Create Omics Object  ###
  obj <- new(
    "OmicsSurv",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls,
    eventTime = eventTime_num,
    eventObserved = eventObserved_lgl
  )

  IntersectOmicsPwyCollct(obj, trim = minPathSize)

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
                            response_num,
                            minPathSize = 3,
                           ...){

  ###  Error Checks and Warnings for supplied assay  ###
  assayData_df <- CheckAssay(assayData_df, ...)


  ###  Process Pathway Collection  ###
  pathwayCollection_ls$setsize <- lengths(pathwayCollection_ls$pathways)


  ###  Pathways List Names Checking  ###
  # If there are no names, create them. If there are missing names, label them.
  if(is.null(names(pathwayCollection_ls$pathways))){

    pathNames <- paste0("path", 1:length(pathwayCollection_ls$pathways))
    names(pathwayCollection_ls$pathways) <- pathNames

  } else if(anyNA(names(pathwayCollection_ls$pathways))){

    missingNm_idx <- which(is.na(names(pathwayCollection_ls$pathways)))
    names(pathwayCollection_ls$pathways)[missingNm_idx] <-
      paste0("noName", missingNm_idx)
    pathNames <- names(pathwayCollection_ls$pathways)

  } else {
    pathNames <- names(pathwayCollection_ls$pathways)
  }

  if(anyDuplicated(names(pathwayCollection_ls$pathways))){

    shortSet_df <- data.frame(
      lapply(pathwayCollection_ls$pathways, "length<-", 1)
    )
    pathNames <- names(shortSet_df)
    names(pathwayCollection_ls$pathways) <- pathNames

  }


  ###  Add Name Key to TERMS and setsize  ###
  names(pathwayCollection_ls$TERMS) <- pathNames
  names(pathwayCollection_ls$setsize) <- pathNames


  ###  Create Omics Object  ###
  obj <- new(
    "OmicsReg",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls,
    response = response_num
  )

  IntersectOmicsPwyCollct(obj, trim = minPathSize)

}




#' @param response_fact A \code{factor} vector of length \eqn{N}: the dependent
#'   variable of a generalized linear regression exercise.
#' @export
#' @rdname CreateOmicsPathway
CreateOmicsCateg <- function(assayData_df,
                              pathwayCollection_ls,
                              response_fact,
                              minPathSize = 3,
                             ...){

  ###  Error Checks and Warnings for supplied assay  ###
  assayData_df <- CheckAssay(assayData_df, ...)


  ###  Process Pathway Collection  ###
  pathwayCollection_ls$setsize <- lengths(pathwayCollection_ls$pathways)


  ###  Pathways List Names Checking  ###
  # If there are no names, create them. If there are missing names, label them.
  #   If there are duplicated names (because R is stupid and allows duplicate
  #   element names in a list -- but not a data frame!), then use the data.frame
  #   name rule to append a period then integers to the end of the name string.
  if(is.null(names(pathwayCollection_ls$pathways))){

    pathNames <- paste0("path", 1:length(pathwayCollection_ls$pathways))
    names(pathwayCollection_ls$pathways) <- pathNames

  } else if(anyNA(names(pathwayCollection_ls$pathways))){

    missingNm_idx <- which(is.na(names(pathwayCollection_ls$pathways)))
    names(pathwayCollection_ls$pathways)[missingNm_idx] <-
      paste0("noName", missingNm_idx)
    pathNames <- names(pathwayCollection_ls$pathways)

  } else {
    pathNames <- names(pathwayCollection_ls$pathways)
  }

  if(anyDuplicated(names(pathwayCollection_ls$pathways))){

    shortSet_df <- data.frame(
      lapply(pathwayCollection_ls$pathways, "length<-", 1)
    )
    pathNames <- names(shortSet_df)
    names(pathwayCollection_ls$pathways) <- pathNames

  }


  ###  Add Name Key to TERMS and setsize  ###
  names(pathwayCollection_ls$TERMS) <- pathNames
  names(pathwayCollection_ls$setsize) <- pathNames


  ###  Create Omics Object  ###
  obj <- new(
    "OmicsCateg",
    assayData_df = assayData_df,
    pathwayCollection = pathwayCollection_ls,
    response = response_fact
  )

  IntersectOmicsPwyCollct(obj, trim = minPathSize)

}
