#' Subset a \code{pathwayCollection}-class Object.
#'
#' @description The subset method for pathways lists as returned by the
#'    \code{\link{read_gmt}} function.
#'
#' @param x An object of class \code{pathwayCollection}.
#' @param name_char The name of a pathway in the collection
#'
#' @return A list of the pathway name (\code{pathway}), the pathway contents
#'    (\code{IDs}), and the pathway description (\code{description}).
#'
#' @details This function finds the index matching the \code{name_char} argument
#'    to the \code{TERMS} field of the \code{pathwayCollection}-class Object,
#'    then subsets the \code{pathways} list, \code{TERMS} vector, and
#'    \code{description} vector by this index. If you subset a trimmed
#'    \code{pathwayCollection} object, and the function errors with "Pathway not
#'    found.", then the pathway specified has been trimmed from the pathway
#'    collection.
#'
#'    Also, this function does not allow for users to overwrite any portion of
#'    a pathway collection. These objects should rarely, if ever, be changed.
#'    If you absolutely must change the components of a \code{pathwayCollection}
#'    object, then create a new one with the code{\link{createPathwayCollection}}
#'    function.
#'
#' @export
#'
#' @examples
#'   data("colon_pathwayCollection")
#'   colon_pathwayCollection["KEGG_RETINOL_METABOLISM"]
#'
`[[.pathwayCollection` <- function(x, name_char){

  path_idx <- which(x$TERMS == name_char)
  if(length(path_idx) == 0) stop("Pathway not found.")

  list(
    Pathway = x$TERMS[path_idx],
    IDs = x$pathways[[path_idx]],
    Description = ifelse(
      test = is.null(x$description),
      yes = NA,
      no = x$description[path_idx]
    )
  )

}
