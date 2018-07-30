#' Display the Summary of a \code{pathwaySet}-class Object.
#'
#' @description The display method for pathways lists as returned by the
#'    \code{\link{read_gmt}} function.
#'
#' @param x An object of class \code{pathwaySet}.
#' @param ... Lazy dots for additional internal arguments (currently unused).
#'
#' @return \code{x}, returned invisibly (with the \code{\link{invisible}}
#'    function).
#'
#' @details This function sets a \code{print} method for \code{pathwaySet}
#'    objects.
#'
#' @export
#'
#' @seealso \code{\link{read_gmt}}; \code{\link{write_gmt}}
#'
#' @importFrom utils str
#'
#' @examples
#'   ###  Load the Example Data  ###
#'   data("colon_pathwaySet")
#'
#'   ###  Print / Show  ###
#'   colon_pathwaySet
#'
print.pathwaySet <- function(x, ...){

  classes_char <- class(x)
  cat("Object with Class(es) '",
      paste(classes_char, collapse = "', '"),
      "' [package 'pathwayPCA'] with ",
      length(x), " elements: \n",
      sep = "")
  str(x,
      max.level = 1,
      vec.len = 1,
      give.attr = FALSE,
      no.list = TRUE)

}
