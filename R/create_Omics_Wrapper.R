#' Generation Wrapper function for \code{-Omics*}-class objects
#'
#' This function calls the \code{\link{create_OmicsPath}},
#'    \code{\link{create_OmicsSurv}}, \code{\link{create_OmicsReg}}, and
#'    \code{\link{create_OmicsCateg}} functions to create valid objects of the
#'    classes \code{OmicsPathway}, \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg}, respectively.
#'
#' @param assayData_df An \eqn{N \times p} data frame with named columns.
#' @param pathwaySet_ls A \code{pathwaySet} list of known gene pathways with two
#'   elements:
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
#'   }
#' @param response An optional response object. See "Details" for more
#'    information. Defaults to \code{NULL}.
#' @param respType What type of response has been supplied. Options are
#'    \code{"none"}, \code{"survival"}, \code{"regression"}, and
#'    \code{"categorical"}. Defaults to \code{"none"} to match the default
#'    \code{response = NULL} value.
#'
#' @details This function is a wrapper around the four \code{create_Omics*}
#'    functions, and as such is not necessary for analysis. This function simply
#'    exists to make \code{Omics*}-object creation easier. The values supplied
#'    to the \code{response} function argument can be in a list, data frame,
#'    matrix, vector, \code{\link[survival]{Surv}} object, or any class which
#'    extends these. Because this function makes "best guess" type conversions
#'    based on the \code{respType} argument, this argument is mandatory if
#'    \code{response} is non-\code{NULL}. Further, it is the responsibility of
#'    the user to ensure that the coerced response contained in the resulting
#'    code{Omics} object accurately reflects the supplied response.
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
#'    \code{\link{create_OmicsPath}},
#'    \code{\link[=OmicsSurv-class]{OmicsSurv}},
#'    \code{\link{create_OmicsSurv}},
#'    \code{\link[=OmicsCateg-class]{OmicsCateg}},
#'    \code{\link{create_OmicsCateg}}
#'    \code{\link[=OmicsReg-class]{OmicsReg}}, and
#'    \code{\link{create_OmicsReg}}.
#'
#' @importFrom survival is.Surv
#'
#' @examples
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwaySet")
#'
#'   ###  Create an OmicsPathway Object  ###
#'   colon_OmicsPath <- create_Omics(assayData_df = colonSurv_df[, -(1:2)],
#'                                   pathwaySet_ls = colon_pathwaySet)
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_Omics(assayData_df = colonSurv_df[, -(1:2)],
#'                                   pathwaySet_ls = colon_pathwaySet,
#'                                   response = colonSurv_df[, 1:2],
#'                                   respType = "surv")
#'
#'   ###  Create an OmicsReg Object  ###
#'   colon_OmicsReg <- create_Omics(assayData_df = colonSurv_df[, -(1:2)],
#'                                  pathwaySet_ls = colon_pathwaySet,
#'                                  response = colonSurv_df$OS_time,
#'                                  respType = "reg")
#'
#'   ###  Create an OmicsCateg Object  ###
#'   colon_OmicsCateg <- create_Omics(assayData_df = colonSurv_df[, -(1:2)],
#'                                    pathwaySet_ls = colon_pathwaySet,
#'                                    response = colonSurv_df$OS_event,
#'                                    respType = "cat")
#' }
#'
#' @export
create_Omics <- function(assayData_df,
                         pathwaySet_ls,
                         response = NULL,
                         respType = c("none", "survival", "regression", "categorical")){

  # browser()

  respType <- match.arg(respType)

  if(respType != "none" && is.null(response)){
    stop(paste0("Response must be specified for type = ", respType))
  }
  if(respType == "none" && !is.null(response)){
    stop("Response type required when a response is given.")
  }

  if(!is.null(response)){
    respClean <- .convertResponse(response, type = respType)
  }

  switch (respType,
          none = {

            message("Creating object of class OmicsPathway.")
            create_OmicsPath(assayData_df = assayData_df,
                             pathwaySet_ls = pathwaySet_ls)

          },
          survival = {

            message("Creating object of class OmicsSurv.")
            create_OmicsSurv(assayData_df = assayData_df,
                             pathwaySet_ls = pathwaySet_ls,
                             eventTime_num = respClean$time,
                             eventObserved_lgl = respClean$dead)

          },
          regression = {

            message("Creating object of class OmicsReg.")
            create_OmicsReg(assayData_df = assayData_df,
                            pathwaySet_ls = pathwaySet_ls,
                            response_num = respClean)

          },
          categorical = {

            message("Creating object of class OmicsCateg.")
            create_OmicsCateg(assayData_df = assayData_df,
                              pathwaySet_ls = pathwaySet_ls,
                              response_fact = respClean)

          }
  )

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
