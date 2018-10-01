#' Access and Edit Response of an \code{OmicsReg} or \code{OmicsReg} Object
#'
#' @description "Get" or "Set" the values of the \code{response_num} or
#'    \code{response_fact} slots of an object of class \code{OmicsReg} or
#'    \code{OmicsReg}, respectively.
#'
#' @param object An object of class \code{\link{OmicsReg-class}} or
#'    \code{\link{OmicsCateg-class}}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return The "get" functions return the objects in the slots specified:
#'    \code{getResponse} returns the \code{response_num} vector from objects of
#'    class \code{OmicsReg} and the \code{response_fact} vector from objects of
#'    class \code{OmicsCateg}. These functions can extract these values from
#'    any valid object of those classes.
#'
#'    The "set" functions enable the user to edit or replace the object in the
#'    \code{response_num} slot for any \code{OmicsReg} object or
#'    \code{response_fact} slot for any \code{OmicsCateg} object, provided that
#'    the new values do not violate the validity check of such an object. See
#'    "Details" for more information.
#'
#' @details These functions can be useful to set or extract the response vector
#'    from an object of class \code{OmicsReg} or \code{OmicsReg}. However, we
#'    recommend that users simply create a new, valid object instead of
#'    modifying an existing one. The validity of edited objects is checked with
#'    their respective \code{\link{ValidOmicsCateg}} or
#'    \code{\link{ValidOmicsReg}} function. Because both classes have a
#'    \code{response} slot, we set this method for the parent class,
#'    \code{\link{OmicsPathway-class}}.
#'
#' @seealso \code{\link{CreateOmicsReg}}, \code{\link{CreateOmicsCateg}}
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
#'   getResponse(colon_OmicsReg)
#'   getResponse(colon_OmicsCateg)
#'
#'   getResponse(colon_OmicsReg) <- newResponse_num
#'   getResponse(colon_OmicsCateg) <- newResponse_fct
#' }



######  Set Generics  ######

#' @export
#' @rdname get_set_OmicsRegCateg
setGeneric("getResponse",
           function(object, ...){
             standardGeneric("getResponse")
           })

#' @export
#' @rdname get_set_OmicsRegCateg
setGeneric("getResponse<-",
           function(object, value){
             standardGeneric("getResponse<-")
           })



###### Set Methods  ######

#' @rdname get_set_OmicsRegCateg
setMethod(f = "getResponse", signature = "OmicsPathway",
          definition = function(object, ...){
            object@response
          })

#' @param value The replacement object to be assigned to the \code{response}
#'    slot.
#'
#' @rdname get_set_OmicsRegCateg
setMethod(f = "getResponse<-", signature = "OmicsPathway",
          definition = function(object, value){

            object@response <- value

            if(validObject(object)){
              return(object)
            }

          })
