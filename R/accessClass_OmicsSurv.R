#' Access and Edit Event Time or Indicator in an \code{OmicsSurv} Object
#'
#' @description "Get" or "Set" the values of the \code{eventTime_num} or
#'    \code{eventObserved_lgl} slots of an object of class \code{OmicsSurv}.
#'
#' @param object An object of class \code{\link{OmicsSurv-class}}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return The "get" functions return the objects in the slots specified:
#'    \code{getEventTime} returns the \code{eventTime_num} vector object and
#'    \code{getEvent} returns the \code{eventObserved_lgl} vector object. These
#'    functions can extract these values from any valid \code{OmicsSurv} object.
#'
#'    The "set" functions enable the user to edit or replace objects in the
#'    \code{eventTime_num} or \code{eventObserved_lgl} slots for any
#'    \code{OmicsSurv} object, provided that the new values do not violate the
#'    validity check of an \code{OmicsSurv} object. See "Details" for more
#'    information.
#'
#' @details These functions can be useful to set or extract the event time or
#'    death indicator from an \code{OmicsSurv} object. However, we recommend
#'    that users simply create a new, valid \code{OmicsSurv} object instead of
#'    modifying an existing one. The validity of edited objects is checked with
#'    the \code{\link{ValidOmicsSurv}} function.
#'
#' @seealso \code{\link{CreateOmics}}
#'
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsSurv.R
#'
#' @importFrom methods setGeneric
#' @importFrom methods setMethod
#' @importFrom methods validObject
#'
#' @examples
#' \dontrun{
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   colon_Omics <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "survival"
#'   )
#'
#'   getEventTime(colon_Omics)
#'   getEvent(colon_Omics)
#'
#'   getEventTime(colon_Omics) <- newTime_num
#'   getEvent(colon_Omics) <- newEvent_lgl
#' }
#'
#' @name SubsetOmicsSurv
#' @rdname get_set_OmicsSurv
NULL



######  Set Generics  ######

#' @export
#' @rdname get_set_OmicsSurv
setGeneric("getEventTime",
           function(object, ...){
             standardGeneric("getEventTime")
           })

#' @export
#' @rdname get_set_OmicsSurv
setGeneric("getEventTime<-",
           function(object, value){
             standardGeneric("getEventTime<-")
           })

#' @export
#' @rdname get_set_OmicsSurv
setGeneric("getEvent",
           function(object, ...){
             standardGeneric("getEvent")
           })

#' @export
#' @rdname get_set_OmicsSurv
setGeneric("getEvent<-",
           function(object, value){
             standardGeneric("getEvent<-")
           })



###### Set Methods  ######

#' @rdname get_set_OmicsSurv
setMethod(f = "getEventTime", signature = "OmicsSurv",
          definition = function(object, ...){
            object@eventTime
          })

#' @param value The replacement object to be assigned to the specified slot.
#'
#' @rdname get_set_OmicsSurv
setMethod(f = "getEventTime<-", signature = "OmicsSurv",
          definition = function(object, value){

            object@eventTime <- value

            if(validObject(object)){
              return(object)
            }

          })

#' @rdname get_set_OmicsSurv
setMethod(f = "getEvent", signature = "OmicsSurv",
          definition = function(object, ...){
            object@eventObserved
          })

#' @rdname get_set_OmicsSurv
setMethod(f = "getEvent<-", signature = "OmicsSurv",
          definition = function(object, value){

            object@eventObserved <- value

            if(validObject(object)){
              return(object)
            }

          })
