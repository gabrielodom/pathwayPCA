#' Generation Wrapper function for \code{-Omics*}-class objects
#'
#' @description This function calls the \code{\link{CreateOmicsPath}},
#'    \code{\link{CreateOmicsSurv}}, \code{\link{CreateOmicsReg}}, and
#'    \code{\link{CreateOmicsCateg}} functions to create valid objects of the
#'    classes \code{OmicsPathway}, \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg}, respectively.
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
#'     \item{\code{TERMS}: }{ A character vector the same length as the
#'        \code{pathways} list with the proper names of the pathways.}
#'     \item{\code{description} : }{An optional character vector the same length
#'        as the \code{pathways} list with additional information about the
#'        pathways.}
#'   }
#'   If your gene pathways list is stored in a \code{.gmt} file, use the
#'   \code{\link{read_gmt}} function to import your pathways list as a
#'   \code{pathwayCollection} list object.
#' @param response An optional response object. See "Details" for more
#'    information. Defaults to \code{NULL}.
#' @param respType What type of response has been supplied. Options are
#'    \code{"none"}, \code{"survival"}, \code{"regression"}, and
#'    \code{"categorical"}. Defaults to \code{"none"} to match the default
#'    \code{response = NULL} value.
#' @param centerScale Should the values in \code{assayData_df} be centered and
#'    scaled? Defaults to \code{TRUE} for centering and scaling, respectively.
#'    See \code{\link{scale}} for more information.
#' @param minPathSize What is the smallest number of genes allowed in each
#'   pathway? Defaults to 3.
#' @param ... Dots for additional arguments passed to the internal
#'   \code{\link{CheckAssay}} function.
#'
#' @details This function is a wrapper around the four \code{CreateOmics*}
#'    functions. The values supplied to the \code{response} function argument
#'    can be in a list, data frame, matrix, vector, \code{\link[survival]{Surv}}
#'    object, or any class which extends these. Because this function makes
#'    "best guess" type conversions based on the \code{respType} argument, this
#'    argument is mandatory if \code{response} is non-\code{NULL}. Further, it
#'    is the responsibility of the user to ensure that the coerced response
#'    contained in the resulting \code{Omics} object accurately reflects the
#'    supplied response.
#'
#'    For \code{respType = "survival"}, \code{response} is assumed to be ordered
#'    by event time, then event indicator. For example, if the response is a
#'    data frame or matrix, this function assumes that the first column is the
#'    time and the second column the death indicator. If the response is a list,
#'    then this function assumes that the first entry in the list is the event
#'    time and the second entry the death indicator. The death indicator must
#'    be a logical or binary (0-1) vector, where 1 or \code{TRUE} represents a
#'    death and 0 or \code{FALSE} represents right-censoring.
#'
#'    Some of the pathways in the supplied pathways list will be removed, or
#'    "trimmed", during object creation. For the pathway-testing methods, these
#'    trimmed pathways will have \eqn{p}-values given as \code{NA}. For an
#'    explanation of pathway trimming, see the documentation for the
#'    \code{\link{IntersectOmicsPwyCollct}} function.
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
#'    \code{\link{CreateOmicsPath}},
#'    \code{\link[=OmicsSurv-class]{OmicsSurv}},
#'    \code{\link{CreateOmicsSurv}},
#'    \code{\link[=OmicsCateg-class]{OmicsCateg}},
#'    \code{\link{CreateOmicsCateg}}
#'    \code{\link[=OmicsReg-class]{OmicsReg}},
#'    \code{\link{CreateOmicsReg}},
#'    \code{\link{CheckAssay}},
#'    \code{\link{CheckPwyColl}}, and
#'    \code{\link{IntersectOmicsPwyCollct}}
#'
#' @importFrom survival is.Surv
#'
#' @examples
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create an OmicsPathway Object  ###
#'   colon_OmicsPath <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection
#'   )
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "surv"
#'   )
#'
#'   ###  Create an OmicsReg Object  ###
#'   colon_OmicsReg <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:2],
#'     respType = "reg"
#'   )
#'
#'   ###  Create an OmicsCateg Object  ###
#'   colon_OmicsCateg <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, c(1,3)],
#'     respType = "cat"
#'   )
#' }
#'
#' @export
CreateOmics <- function(assayData_df,
                        pathwayCollection_ls,
                        response = NULL,
                        respType = c("none", "survival", "regression", "categorical"),
                        centerScale = c(TRUE, TRUE),
                        minPathSize = 3,
                        ...){

  # browser()

  ###  Match and Check Response  ###
  respType <- match.arg(respType)
  if(respType != "none" && is.null(response)){
    stop(paste0("Response must be specified for type = ", respType))
  }
  if(respType == "none" && !is.null(response)){
    stop("Response type required when a response is given.")
  }


  ###  Data Error Checks and Warnings  ###
  # Assay
  origClass    <- class(assayData_df)
  assayData_df <- CheckAssay(assayData_df, ...)
  assayData_df <- CheckSampleIDs(assayData_df)

  # Pathway Collection
  pathwayCollection_ls <- CheckPwyColl(pathwayCollection_ls)

  # merge the pheno and assay data
  if(!is.null(response)){

    respClean_df <- .convertPhenoDF(response, type = respType)
    respClean_df <- CheckSampleIDs(respClean_df)
    data_ls <- JoinPhenoAssay(
      pheno_df = respClean_df,
      assay_df = assayData_df
    )

    assayData_df <- data_ls$assay
    respClean_df <- data_ls$response
    sampleID     <- data_ls$sampleID

    respClean <- .convertResponse(respClean_df, type = respType)

  } else {

    sampleID     <- assayData_df[, 1, drop = TRUE]
    assayData_df <- assayData_df[, -1]

  }


  ###  Centre and Scale Assay  ###
  if(any(centerScale)){

    assayData_df <- as.data.frame(
      scale(assayData_df, center = centerScale[1], scale = centerScale[2])
    )
    class(assayData_df) <- origClass

  }


  ###  Create Data Object  ###
  switch (respType,
          none = {

            message("\n  ======  Creating object of class OmicsPathway  ======")
            out <- CreateOmicsPath(
              assayData_df = assayData_df,
              sampleIDs_char = sampleID,
              pathwayCollection_ls = pathwayCollection_ls
            )

          },
          survival = {

            message("\n  ======  Creating object of class OmicsSurv  =======")
            out <- CreateOmicsSurv(
              assayData_df = assayData_df,
              sampleIDs_char = sampleID,
              pathwayCollection_ls = pathwayCollection_ls,
              eventTime_num = respClean$time,
              eventObserved_lgl = respClean$dead
            )

          },
          regression = {

            message("\n  ======  Creating object of class OmicsReg  =======")
            out <- CreateOmicsReg(
              assayData_df = assayData_df,
              sampleIDs_char = sampleID,
              pathwayCollection_ls = pathwayCollection_ls,
              response_num = respClean
            )

          },
          categorical = {

            message("\n  ======  Creating object of class OmicsCateg  =======")
            out <- CreateOmicsCateg(
              assayData_df = assayData_df,
              sampleIDs_char = sampleID,
              pathwayCollection_ls = pathwayCollection_ls,
              response_fact = respClean
            )

          }
  )

  ###  Return  ###
  IntersectOmicsPwyCollct(out, trim = minPathSize)

}

.convertResponse <- function(object,
                             type = c("survival", "regression", "categorical")){

  type <- match.arg(type)

  if(type == "survival"){

    if(is.Surv(object)){
      object <- as.matrix(object)
    }

    if(is.data.frame(object)){
      object <- as.matrix(object)
    }

    if(is.matrix(object)){

      if(ncol(object) == 2){

        outTime <- object[, 1, drop = TRUE]
        outDeath <- object[, 2, drop = TRUE]

      } else {
        stop("Object must have two columns only: death time and death indicator.")
      }

    } else if(is.list(object)){

      if(length(object) == 2){

        outTime <- object[[1]]
        outDeath <- object[[2]]

      } else {
        stop("Object must have two entries only: death time and death indicator.")
      }
    }

    if(is.atomic(outTime) && is.atomic(outDeath)){

      outTime <- as.numeric(outTime)
      outDeath <- as.logical(outDeath)

    } else {
      stop("Death time and death indicator must be atomic vectors.")
    }

    list(time = outTime,
         dead = outDeath)

  } else {

    if(is.data.frame(object)){
      object <- as.matrix(object)
    }

    if(is.matrix(object)){

      if(ncol(object) == 1 || nrow(object) == 1){
        object <- as.vector(object)
      } else {
        stop("Multivariate response not supported at this time.")
      }

    }

    if(is.list(object)){

      if(length(object) == 1){
        object <- unlist(object, use.names = FALSE)
      } else {
        stop("Non-survival response must be a vector or 1-dimensional data frame.")
      }
    }

    if(is.atomic(object)){

      if(type == "regression"){
        object <- as.numeric(object)
      } else {
        object <- as.factor(object)
      }

    } else {
      stop("Non-survival response must be a vector or 1-dimensional data frame.")
    }

    object

  }

}


.convertPhenoDF <- function(pheno_df,
                            type = c("survival", "regression", "categorical")){

  type <- match.arg(type)

  if(!inherits(pheno_df, "data.frame")){
    stop("The response must be a data frame.")
  }

  if(type == "survival"){

    if(ncol(pheno_df) != 3){
      stop("Survival data must be a data frame with three columns, sample ID, event time,
  and death indicator, in exactly that order.")
    } else {

      pheno_df[, 2] <- as.numeric(pheno_df[, 2, drop = TRUE])
      pheno_df[, 3] <- as.logical(pheno_df[, 3, drop = TRUE])

    }

  } else {

    if(ncol(pheno_df) != 2){
      stop("Regression and categorical data must be a data frame with two columns, sample ID
  and response, in exactly that order.")
    } else {

      if(type == "regression"){
        pheno_df[, 2] <- as.numeric(pheno_df[, 2, drop = TRUE])
      } else {
        pheno_df[, 2] <- as.factor(pheno_df[, 2, drop = TRUE])
      }

    }

  }

  pheno_df

}
