#' Access and Edit Assay or \code{pathwayCollection} Values in \code{Omics*}
#'    Objects
#'
#' @description "Get" or "Set" the values of the \code{assayData_df},
#'    \code{sampleIDs_char}, or \code{pathwayCollection} slots of an object of
#'    class \code{OmicsPathway} or a class that extends this class
#'    (\code{OmicsSurv}, \code{OmicsReg}, or \code{OmicsCateg}).
#'
#' @param object An object of or extending \code{\link{OmicsPathway-class}}:
#'    that class, \code{\link{OmicsSurv-class}}, \code{\link{OmicsReg-class}},
#'    or \code{\link{OmicsCateg-class}}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return The "get" functions return the objects in the slots specified:
#'    \code{getAssay} returns the \code{assayData_df} data frame object,
#'    \code{getSampleIDs} returns the \code{sampleIDs_char} character vector,
#'    \code{getPathwayCollection} returns the \code{pathwayCollection} list
#'    object, and \code{getTrimPathwayCollection} returns the
#'    \code{trimPathwayCollection}. These functions can extract these values
#'    from any valid \code{OmicsPathway}, \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg} object.
#'
#'    The "set" functions enable the user to edit or replace objects in the
#'    \code{assayData_df}, \code{sampleIDs_char}, or \code{pathwayCollection}
#'    slots for any \code{OmicsPathway}, \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg} objects, provided that the new values do not violate
#'    the validity checks of their respective objects. Because the slot for
#'    \code{trimPathwayCollection} is filled upon object creation, and to ensure
#'    that this pathway collection is "clean", there is no "set" function for
#'    the trimmed pathway collection slot. Instead, users can update the pathway
#'    collection, and the trimmed pathway collection will be updated
#'    automatically. See "Details" for more information on the "set" functions.
#'
#' @details These functions can be useful to set or extract the assay data or
#'    pathways list from an \code{Omics*}-class object. However, we recommend
#'    that users simply create a new, valid \code{Omics*} object instead of
#'    modifying an existing one. The validity of edited objects is checked with
#'    the \code{\link{ValidOmicsSurv}}, \code{\link{ValidOmicsCateg}}, or
#'    \code{\link{ValidOmicsReg}} functions.
#'
#'    Further, because the \code{pathwayPCA} methods require a cleaned (trimmed)
#'    pathway collection, the \code{trimPathwayCollection} slot is read-only.
#'    Users may only edit this slot by updating the pathway collection provided
#'    to the \code{pathwayCollection} slot. Despite this functionality, we
#'    \strong{strongly} recommend that users create a new object with the
#'    updated pathway collection, rather than attempting to overwrite the slots
#'    within an existing object. See \code{\link{IntersectOmicsPwyCollct}} for
#'    details on trimmed pathway collection.
#'
#' @seealso \code{\link{CreateOmics}}
#'
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @importFrom methods setGeneric
#' @importFrom methods setMethod
#' @importFrom methods validObject
#'
#'
#' @examples
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   colon_Omics <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection
#'   )
#'
#'   getAssay(colon_Omics)
#'   getPathwayCollection(colon_Omics)
#'
#' @name SubsetOmicsPath
#' @rdname get_set_OmicsPathway
NULL



######  Create Generics  ######

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
setGeneric("getSampleIDs",
           function(object, ...){
             standardGeneric("getSampleIDs")
           })

#' @export
#' @rdname get_set_OmicsPathway
setGeneric("getSampleIDs<-",
           function(object, value){
             standardGeneric("getSampleIDs<-")
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

#' @export
#' @rdname get_set_OmicsPathway
setGeneric("getTrimPathwayCollection",
           function(object, ...){
             standardGeneric("getTrimPathwayCollection")
           })



######  Create Methods  ######

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
setMethod(f = "getSampleIDs", signature = "OmicsPathway",
          definition = function(object, ...){
            object@sampleIDs_char
          })

#' @rdname get_set_OmicsPathway
setMethod(f = "getSampleIDs<-", signature = "OmicsPathway",
          definition = function(object, value){

            object@sampleIDs_char <- value

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
              return(IntersectOmicsPwyCollct(object))
            }

          })

#' @rdname get_set_OmicsPathway
setMethod(f = "getTrimPathwayCollection", signature = "OmicsPathway",
          definition = function(object, ...){
            object@trimPathwayCollection
          })
