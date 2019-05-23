#' Filter and Subset a \code{pathwayCollection}-class Object by Symbol.
#'
#' @description The filter-subset method for pathways lists as returned by the
#'    \code{\link{read_gmt}} function. This function returns the subset of 
#'    pathways which contain the set of symbols requested
#'
#' @param x An object of class \code{pathwayCollection}.
#' @param symbols_char A character vector or scalar of gene symbols or regions
#' @param ... Additional arguments passed to the \code{\link{Contains}} function
#'
#' @return An object of class \code{pathwayCollection}, but containing only the
#'    sets which contain the symbols supplied to \code{symbols_char}. If no sets
#'    are found to contain the symbols supplied, this function returns
#'    \code{NULL} and prints a warning.
#'
#' @details This function finds the index of each set that contains the symbols
#'    supplied, then returns those sets as a new \code{pathwayCollection}
#'    object. Find pathways that contain geneA OR geneB by passing the argument
#'    \code{matches = "any"} through \code{...} to \code{\link{Contains}} (this
#'    is the default value). Find pathways that contain geneA AND geneB by
#'    changing this argument to \code{matches = "all"}. Find all genes in a 
#'    specified family by passing in one value to \code{short} and setting
#'    \code{partial = TRUE}.
#'
#' @export
#'
#' @examples
#'   data("colon_pathwayCollection")
#'   
#'   WhichPathways(colon_pathwayCollection, "MAP", partial = TRUE)
#'   
#'   WhichPathways(
#'     colon_pathwayCollection,
#'     c("MAP4K5", "RELA"),
#'     matches = "all"
#'   )
#'
WhichPathways <- function(x, symbols_char, ...){
  
  if(!("pathwayCollection" %in% class(x))){
    stop("Object must be or extend class 'pathwayCollection'.", call. = FALSE)
  }
  type_char <- names(x)[1]
  
  symbol_logi <- sapply(
    x[[1]],
    Contains,
    short = symbols_char,
    ...
  )
  
  if(sum(symbol_logi) == 0){
    warning("No sets found to contain the requested symbols.", call. = FALSE)
    return(NULL)
  }
  
  if(is.null(x$description)){
    CreatePathwayCollection(
      sets_ls = x[[1]][symbol_logi],
      TERMS = x$TERMS[symbol_logi],
      setType = type_char
    )
  } else {
    CreatePathwayCollection(
      sets_ls = x[[1]][symbol_logi],
      TERMS = x$TERMS[symbol_logi],
      description = x$description[symbol_logi],
      setType = type_char
    )
  }
  
}
