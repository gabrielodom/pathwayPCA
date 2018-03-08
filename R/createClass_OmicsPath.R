#' An S4 class for mass spectronomy data and gene pathway lists
#'
#' @slot massSpec An $N x p$ data frame with named columns
#' @slot pathwaySet A list of known gene pathways with one or two elements:
#' \itemize{
#'   \item{pathways : }{A named list of character vectors. Each vector contains
#'     the names of the individual genes within that pathway as a vector of
#'     character strings. The names contained in these vectors must have non-
#'     empty overlap with the \emph{column names} of the \code{massSpec} data
#'     frame. The names of the pathways (the list elements themselves) should
#'     be the a shorthand representation of the full pathway name.}
#'   \item{TERMS: }{ A character vector the same length as the
#'     \code{pathways} list with the proper names of the pathways.}
#'   \item{setsize : }{A named integer vector the same length as the
#'     \code{pathways} list with the number of genes in each pathway. This list
#'     item is calculated during the creation step of a \code{create_Omics*()}
#'     function call.}
#'
#' }
#'
#' @seealso \code{\link{create_OmicsPath}}
#'
#' @importFrom methods new
#'
#' @export
setClass("OmicsPathway",
         slots = c(massSpec = "data.frame", pathwaySet = "list"))
# so my create.*() functions don't throw a fit when they see a tibble:
setOldClass(c("tbl_df", "tbl", "data.frame"))
