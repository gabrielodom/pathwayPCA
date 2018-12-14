#' Subset a \code{pathwayCollection}-class Object.
#'
#' @description The subset method for pathways lists as returned by the
#'    \code{\link{read_gmt}} function.
#'
#' @param x An object of class \code{pathwayCollection}.
#' @param name_char The name of a pathway in the collection or its unique ID.
#'
#' @return A list of the pathway name (\code{Term}), unique ID (\code{pathID}),
#'    contents (\code{IDs}), description (\code{description}), and number of
#'    features (\code{Size}).
#'
#' @details This function finds the index matching the \code{name_char} argument
#'    to the \code{TERMS} field of the \code{pathwayCollection}-class Object,
#'    then subsets the \code{pathways} list, \code{TERMS} vector,
#'    \code{description} vector, and \code{setsize} vector by this index. If you
#'    subset a trimmed \code{pathwayCollection} object, and the function errors
#'    with "Pathway not found.", then the pathway specified has been trimmed
#'    from the pathway collection.
#'
#'    Also, this function does not allow for users to overwrite any portion of
#'    a pathway collection. These objects should rarely, if ever, be changed.
#'    If you absolutely must change the components of a \code{pathwayCollection}
#'    object, then create a new one with the code{\link{CreatePathwayCollection}}
#'    function.
#'
#' @export
#'
#' @examples
#'   data("colon_pathwayCollection")
#'   colon_pathwayCollection[["KEGG_RETINOL_METABOLISM"]]
#'
#' @name SubsetPathwayCollection
`[[.pathwayCollection` <- function(x, name_char){
  # browser()

  ###  Checks  ###
  y <- NextMethod("[[")
  if(missing(name_char)) return(y)
  if(!is.character(name_char)) return(y)


  ###  Name Match  ###
  term_idx <- which(x$TERMS == name_char)
  pathID_idx <- which(names(x$pathways) == name_char)
  if(length(term_idx) > 0){
    path_idx <- term_idx
  } else if(length(pathID_idx) > 0){
    path_idx <- pathID_idx
  } else {
    stop("Pathway not found.")
  }


  ###  Calculations and Subsettings  ###
  features_char <- x$pathways[[path_idx]]
  desc_char <- ifelse(is.null(x$description), NA, x$description[path_idx])


  ###  Return  ###
  list(
    Term = x$TERMS[path_idx],
    PathID = names(x$pathways)[path_idx],
    IDs = features_char,
    Description = desc_char,
    Size = length(features_char)
  )

}
