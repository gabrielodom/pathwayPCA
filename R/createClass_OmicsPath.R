#' An S4 class for mass spectrometry or bio-assay data and gene pathway lists
#'
#' @slot assayData_df An \eqn{N \times p} data frame with named columns.
#' @slot pathwayCollection A list of known gene pathways with three or four
#'    elements:
#' \itemize{
#'   \item{\code{pathways} : }{A named list of character vectors. Each vector
#'      contains the names of the individual genes within that pathway as a
#'      vector of character strings. The names contained in these vectors must
#'      have non-empty overlap with the \emph{column names} of the
#'      \code{assayData_df} data frame. The names of the pathways (the list
#'      elements themselves) should be the a shorthand representation of the
#'      full pathway name.}
#'   \item{\code{TERMS} : }{A character vector the same length as the
#'     \code{pathways} list with the proper names of the pathways.}
#'   \item{\code{description} : }{An optional character vector the same length
#'      as the \code{pathways} list with additional information about the
#'      pathways.}
#'   \item{\code{setsize} : }{A named integer vector the same length as the
#'     \code{pathways} list with the number of genes in each pathway. This list
#'     item is calculated during the creation step of a \code{create_OmicsReg}
#'     function call.}
#' }
#' @slot trimPathwayCollection A subset of the list stored in the
#'    \code{pathwayCollection} slot. This list will have pathways that only
#'    contain genes that are present in the assay data frame.
#'
#' @seealso \code{\link{create_OmicsPath}}
#'
#' @importFrom methods setClass
#' @importFrom methods setOldClass
#'
#' @export
setClass(
  "OmicsPathway",
  slots = c(assayData_df = "data.frame",
            pathwayCollection = "list",
            trimPathwayCollection = "list")
)
# so my create.*() functions don't throw a fit when they see a tibble:
setOldClass(c("tbl_df", "tbl", "data.frame"))
# So the lists imported with the read_gmt function are allowed
setOldClass(c("pathwayCollection", "list"))
