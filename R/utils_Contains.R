#' Check if a long atomic vector contains a short atomic vector
#'
#' @description Check if any or all of the elements of a short atomic vector
#'    are contained within a supplied long atomic vector.
#'
#' @param long A vector to possibly containing any or all elements of
#'    \code{short}
#' @param short A short vector or scalar, some elements of which may be
#'    contained in \code{long}
#' @param matches Should partial set matching of \code{short} be allowed?
#'    Defaults to \code{"any"}, signifying that the function should return
#'    \code{TRUE} if any of the elements of \code{short} are contained in
#'    \code{long}. The other option is \code{"all"}.
#' @param partial Should partial string matching be allowed? Defaults to
#'    \code{FALSE}. Partial string matching means that the character string
#'    \strong{starts with} the supplied value.
#'
#' @details This is a helper function to find out if a gene symbol or some
#'    similar character string (or character vector) is contained in a pathway.
#'    Currently, this function uses base R, but we can write it in a compiled
#'    language (such as C++) to increase speed later.
#'    
#'    For partial matching (\code{partial = TRUE}), \code{long} must be an
#'    atomic vector of type character, \code{short} must be an atomic scalar (a
#'    vector with length of 1) of type character, and \code{matches} should be
#'    set to \code{"any"}. Because this function is designed to match gene
#'    symbols or CpG locations, we care if the symbol or location starts with
#'    the string supplied. For example, if we set \code{short = "PIK"}, then we
#'    want to find if any of the gene symbols in the supplied \code{long}
#'    vector belong to the PIK gene family. We don't care if this string appears
#'    elsewhere in a gene symbol.
#'
#' @return A logical scalar. If \code{matches = "any"}, this indicates if any
#'    of the elements of \code{short} are contained in \code{long}. If
#'    \code{matches = "all"}, this indicates if all of the elements of
#'    \code{short} are contained in \code{long}. If \code{partial = TRUE}, the
#'    returned logical indicates whether or not any of the character strings in
#'    \code{long} \strong{start with} the character scalar supplied to
#'    \code{short}.
#'
#' @export
#'
#' @examples
#'    Contains(1:10, 8)
#'    Contains(LETTERS, c("A", "!"), matches = "any")
#'    Contains(LETTERS, c("A", "!"), matches = "all")
#'    
#'    genesPI <- c(
#'      "PI4K2A", "PI4K2B", "PI4KA", "PI4KB", "PIK3C2A", "PIK3C2B", "PIK3C2G",
#'      "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2",
#'      "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PIKFYVE", "PIP4K2A",
#'      "PIP4K2B", "PIP5K1B", "PIP5K1C", "PITPNB"
#'    )
#'    Contains(genesPI, "PIK3", partial = TRUE)
#'
Contains <- function(long, short, matches = c("any", "all"), partial = FALSE){
  
  matches <- match.arg(matches)

  if(!partial){
    
    switch(matches,
           any = {
             any(match(long, short, nomatch = 0))
           },
           all = {
             sum(short %in% long) == length(short)
           }
    )
    
  } else {
    
    if(!is.character(short) || !is.character(long)){
      stop(
        "Partial matching permitted only for character vectors.",
        call. = FALSE
      )
    }
    
    if(length(short) > 1){
      stop(
        "Partial matching permitted only for single values of 'short'.",
        call. = FALSE
      )
    }
    
    if(matches == "all"){
      warning(
        "Partial matching takes any matches. matches = 'all' will be ignored.",
        call. = FALSE
      )
    }
    
    any(startsWith(long, short))
    
  }
  
}
