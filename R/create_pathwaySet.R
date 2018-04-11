#' Manually Create a \code{pathwaySet}-class Object.
#'
#' @description Manually create a \code{pathwaySet} list similar to the output
#'    of the \code{\link{read_gmt}} function.
#'
#' @param pathways A named list of character vectors. Each vector should contain
#'    the names of the individual genes or proteins within that pathway as a
#'    vector of character strings. The names contained in these vectors should
#'    have non-empty overlap with the feature names of the assay data frame that
#'    will be paired with this list in the subsequent analysis. The names of the
#'    pathways (the list elements themselves) should be the a shorthand
#'    representation of the full pathway name, but this is not required.
#' @param TERMS A character vector the same length as the \code{pathways} list
#'    with the proper names of the pathways.
#' @param ... Additional vectors or data components related to the
#'    \code{pathways} list. These values should be passed as a name-value pair.
#'    See "Details" for more information.
#'
#' @return A list object with class \code{pathwaySet}.
#'
#' @details This function checks the pathwy list and pathway term inputs and
#'    then creates a \code{pathwaySet} object from them. Pass additional list
#'    elements using the form \code{tag = value} through the \code{...}
#'    argument (as in the \code{\link{list}} function). Because some functions
#'    in the \code{pathwayPCA} package add and edit elements of
#'    \code{pathwaySet} objects, please do not create \code{pathwaySet} list
#'    items named \code{setsize} or \code{trim_setsize}.
#'
#' @export
#'
#' @seealso \code{\link{read_gmt}}
#'
#' @examples
#'   data("colon_pathwaySet")
#'
#'   create_pathwaySet(pathways = colon_pathwaySet$pathways,
#'                     TERMS = colon_pathwaySet$TERMS)
#'
create_pathwaySet <- function(pathways, TERMS, ...){

  ###  Class Checks  ###
  if(!is.list(pathways)){
    stop("The pathways object must be a list of pathway vectors.")
  }

  if(!is.atomic(TERMS)){
    stop("The TERMS object must be an atomic vector.")
  }

  ###  Validation Checks  ###
  nPaths <- length(pathways)
  if(length(TERMS) != nPaths){
    stop("The TERMS vector must have the same length as the pathways list.")
  }

  dotNames <- names(list(...))
  if(any(c("setsize", "trim_setsize") %in% dotNames)){
    warning("The names 'setsize' and 'trim_setsize' are reserved names of a pathwaySet object.
            Values stored with these names may be overwritten during pathwayPCA function execution.
            Use with extreme caution.")
  }

  ###  Create and Return pathwaySet  ###
  out_ls <- list(pathways = pathways,
                 TERMS = TERMS,
                 ...)
  class(out_ls) <- c("pathwaySet", "list")

  out_ls

  }
