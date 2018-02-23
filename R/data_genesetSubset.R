#'  Gene Pathway Sets
#'
#' @description A Canonical Pathways Gene Sets from the Broad Institute: File:
#'   c2.cp.v6.0.symbols.gmt.
#'
#' @details This is a subset of 15 pathways from the Broad Institute gene set
#'   list. This subset contains seven pathways which are related to the response
#'   information in the \code{\link{colonSurv_df}}
#'
#' @format A list of three elements:
#' \itemize{
#'   \item{pathways : }{A list of 15 character vectors. Each vector contains
#'     the names of the individual genes within that pathway as a vector of
#'     character strings.}
#'   \item{TERMS : }{A character vector of length 15 containing the names of
#'     the gene pathways.}
#'   \item{setsize : }{A named integer vector of length 15 containing the
#'     number of genes in each gene pathway.}
#' }
#'
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
"colonGenesets_ls"