#'  Gene Pathway Sets for Supervised Example
#'
#' @description A Gene Ontology Pathways Gene Set. This is a list of genes with
#'    Entrez IDs and GO annotations (\code{hgu133plus2ENTREZID} and
#'    \code{hgu133plus2GO} from \code{Bioconductor} package
#'    \code{hgu133plus2.db}, respectively), with pathways with only one gene or
#'    more than 180 genes removed (99.5\% of all pathways we encounter have less
#'    than 180 genes).
#'
#' @format A list of three elements:
#' \itemize{
#'   \item{pathways : }{A list of 7949 character vectors. Each vector contains
#'     the ID numbers (as a character) of the individual genes within that
#'     pathway as a vector of character strings.}
#'   \item{TERMS : }{A character vector of length 7949 containing the names of
#'     the gene pathways.}
#'   \item{setsize : }{A named integer vector of length 7949 containing the
#'     number of genes in each gene pathway.}
#' }
#'
#' @source UNKNOWN - GET INFO FROM STEVEN
"supervised_Genesets_ls"
