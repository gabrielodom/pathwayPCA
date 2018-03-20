#'  Gene Pathway Subset
#'
#' @description An example Canonical Pathways Gene Subset from the Broad
#'   Institute: File: \code{c2.cp.v6.0.symbols.gmt}.
#'
#' @details This is a subset of 15 pathways from the Broad Institute gene set
#'   list. This subset contains seven pathways which are related to the response
#'   information in the \code{\link{colonSurv_df}} data file.
#'
#' @format A list of two elements:
#' \itemize{
#'   \item{\code{pathways} : }{A list of 15 character vectors. Each vector
#'      contains the names of the individual genes within that pathway as a
#'      vector of character strings.}
#'   \item{\code{TERMS} : }{A character vector of length 15 containing the
#'      names of the gene pathways.}
#' }
#'
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
"colonGenesets_ls"
