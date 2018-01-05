#'  Raw and Adjusted p-Values for Gene Pathway Sets from Supervised PCA
#'
#' @description The raw and adjusted p-values corresponding to the association
#'    between the gene pathways in the \code{\link{supervised_Tumors_df}} data
#'    set and the survival outcomes in the \code{\link{supervised_patInfo_df}}
#'    data set. The gene pathways have GO annotations (from \code{hgu133plus2GO}
#'    from the \code{Bioconductor} package \code{hgu133plus2.db}). Adjusted
#'    p-values were calculated via the Benjamini & Hochberg (1995) step-up False
#'    Discovery Rate (FDR)-controlling procedure in the
#'    \code{\link{adjustRaw_pVals}} function.
#'
#' @format A data frame with 4,240 rows and five columns:
#' \itemize{
#'   \item{goterms : }{The GO identifiers for each of the 4,240 pathways
#'     considered.}
#'   \item{setsize : }{The number of genes in each gene pathway.}
#'   \item{rawp : }{The pathway p-values before adjustement for FDR.
#'     ???? MORE INFO NEEDED.}
#'   \item{FDR : }{The p-values for each pathway after adjustment to account
#'     for FDR.}
#'   \item{terms : }{A short description of the biological functionality of each
#'     pathway.}
#' }
#'    The data frame's rows are sorted by FDR, with ties broken by the raw
#'    p-values.
#'
#' @source Calculated via the \code{\link{pathway_pValues}} function, which is
#'    a wrapper function for \code{\link{adjustRaw_pVals}}.
"spcaPathwayPvals_df"
