#' Access and Edit Assay or \code{pathwayCollection} Values in \code{Omics*}
#'    Objects
#'
#' @description "Get" or "Set" the values of the \code{assayData_df} or
#'    \code{pathwayCollection} slots of an object of class \code{OmicsPathway}
#'    or a class that extends this class (\code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg}).
#'
#' @param object An object of or extending \code{\link{OmicsPathway-class}}:
#'    that class, \code{\link{OmicsSurv-class}}, \code{\link{OmicsReg-class}},
#'    or \code{\link{OmicsCateg-class}}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return The "get" functions return the objects in the slots specified:
#'    \code{getAssay} returns the \code{assayData_df} data frame object and
#'    \code{getPathwayCollection} returns the \code{pathwayCollection} list
#'    object. These functions can extract these values from any valid
#'    \code{OmicsPathway}, \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg} object.
#'
#'    The "set" functions enable the user to edit or replace objects in the
#'    \code{assayData_df} or \code{pathwayCollection} slots for any
#'    \code{OmicsPathway}, \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg} objects, provided that the new values do not violate
#'    the validity checks of their respective objects. See "Details" for more
#'    information.
#'
#' @details These functions can be useful to set or extract the assay data or
#'    pathways list from an \code{Omics*}-class object. However, we recommend
#'    that users simply create a new, valid \code{Omics*} object instead of
#'    modifying an existing one. The validity of edited objects is checked with
#'    the \code{\link{valid_OmicsSurv}}, \code{\link{valid_OmicsCateg}}, or
#'    \code{\link{valid_OmicsReg}} functions.
#'
#' @seealso \code{\link{create_OmicsPath}}, \code{\link{create_OmicsSurv}},
#'    \code{\link{create_OmicsReg}}, \code{\link{create_OmicsCateg}}
#'
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @importFrom methods setGeneric
#' @importFrom methods setMethod
#' @importFrom methods validObject
#'
#' @examples
#' \dontrun{
#'   getAssay(colon_OmicsSurv)
#'   getPathwayCollection(colon_OmicsSurv)
#'
#'   getAssay(colon_OmicsSurv) <- newAssay_df
#'   getPathwayCollection(colon_OmicsSurv) <- new_pathwayCollection
#' }



######  Set Generics  ######

#' @export
#' @rdname get_set_OmicsPathway
setGeneric("getAssay",
           function(object, ...){
             standardGeneric("getAssay")
           })

#' @export
#' @rdname get_set_OmicsPathway
setGeneric("getAssay<-",
           function(object, value){
             standardGeneric("getAssay<-")
           })

#' @export
#' @rdname get_set_OmicsPathway
setGeneric("getPathwayCollection",
           function(object, ...){
             standardGeneric("getPathwayCollection")
           })

#' @export
#' @rdname get_set_OmicsPathway
setGeneric("getPathwayCollection<-",
           function(object, value){
             standardGeneric("getPathwayCollection<-")
           })



######  Set Methods  ######

#' @rdname get_set_OmicsPathway
setMethod(f = "getAssay", signature = "OmicsPathway",
          definition = function(object, ...){
            object@assayData_df
          })


#' @param value The replacement object to be assigned to the specified slot.
#'
#' @rdname get_set_OmicsPathway
setMethod(f = "getAssay<-", signature = "OmicsPathway",
          definition = function(object, value){

            object@assayData_df <- value

            if(validObject(object)){
              return(object)
            }

          })

#' @rdname get_set_OmicsPathway
setMethod(f = "getPathwayCollection", signature = "OmicsPathway",
          definition = function(object, ...){
            object@pathwayCollection
          })

#' @rdname get_set_OmicsPathway
setMethod(f = "getPathwayCollection<-", signature = "OmicsPathway",
          definition = function(object, value){

            object@pathwayCollection <- value

            if(validObject(object)){
              return(object)
            }

          })
