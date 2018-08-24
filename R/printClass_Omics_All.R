#' Display the Summary of an \code{Omics*}-class Object.
#'
#' @description The display method for objects of class \code{OmicsPathway},
#'    \code{OmicsSurv}, \code{OmicsReg}, or \code{OmicsCateg}.
#'
#' @param object An object inheriting the super-class \code{OmicsPathway}. This
#'    class includes objects of class \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg}.
#'
#' @return A copy of \code{object}, returned invisibly (with the
#'    \code{\link{invisible}} function).
#'
#' @details S4 objects print to the screen via the \code{\link[methods]{show}}
#'    function. This function sets a \code{show} method for \code{OmicsPathway}
#'    objects.
#'
#' @export
#'
#' @include createClass_OmicsPath.R
#'
#' @importFrom methods setMethod
#' @importFrom methods show
#' @importFrom utils str
#'
#' @examples
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_OmicsSurv(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     eventTime_num = colonSurv_df$OS_time,
#'     eventObserved_lgl = as.logical(colonSurv_df$OS_event)
#'   )
#'
#'   ###  Print / Show  ###
#'   colon_OmicsSurv
#'
setMethod(f = "show", signature = "OmicsPathway",
          definition = function(object){
            str(object, max.level = 2)
          })
