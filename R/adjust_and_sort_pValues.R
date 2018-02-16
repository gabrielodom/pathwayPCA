#' Adjust and Sort Pathway \eqn{p}-Values
#'
#' @description Adjust the pathway \eqn{p}-values, then return a data frame of
#'   the relevant pathway information, sorted by adjusted significance.
#'
#' @param pVals_vec A named vector of permutation \eqn{p}-values returned by the
#'   \code{\link{permTest_OmicsSurv}}, \code{\link{permTest_OmicsReg}}, or
#'   \code{\link{permTest_OmicsCateg}} functions when the analysis performed was
#'   AES-PCA. Otherwise, when the analysis was performed with Supervised PCA, a
#'   named vector of \eqn{p}-values from the \code{\link{weibullMix_pValues}}
#'   function.
#' @param genesets_ls A list of known gene pathways. This pathway list must
#'   contain
#'   \itemize{
#'     \item{pathways : }{A named list of character vectors where each vector
#'       contains the names of the genes in that specific pathway.}
#'     \item{TERMS : }{A character vector the same length as \code{pathways}
#'       containing the full pathway descriptions.}
#'     \item{setsize : }{An integer vector the same length as \code{pathways}
#'       containing the number of genes present in the pathway.}
#'   }
#' @param adjust Should you adjust the \eqn{p}-values for multiple comparisons?
#'   Defaults to TRUE.
#' @param proc_vec Character vector of procedures. The returned data frame will
#'   be sorted in ascending order by the first procedure in this vector, with
#'   ties broken by the unadjusted \eqn{p}-value. If only one procedure is
#'   selected, then it is necessarily the first procedure.
#' @param ... Additional arguments to pass to the \code{\link{adjustRaw_pVals}}
#'   function.
#'
#' @return A data frame with columns
#' \itemize{
#'   \item{\code{pathways} : }{The names of the pathways in the \code{Omics*}}
#'     object (stored in \code{object@@pathwaySet$pathways})
#'   \item{\code{setsize} : }{The number of genes in each of the original
#'     pathways (as stored in the \code{object@@pathwaySet$setsize} object)}
#'   \item{\code{terms} : }{The pathway description, as stored in the
#'     \code{object@@pathwaySet$TERMS} object}
#'   \item{\code{rawp} : }{The unadjusted \eqn{p}-values of each pathway}
#'   \item{\code{...} : }{Additional columns as specified through the
#'     \code{adjustment} argument}
#' }
#'
#' The data frame will be sorted in ascending order by the method specified
#'   first in the \code{adjustment} argument. If \code{adjustpValues = FALSE},
#'   then the data frame will be sorted by the raw \eqn{p}-values.
#'
#' @details This is a wrapper function for the \code{\link{adjustRaw_pVals}}
#'   function.
#'
#' @export
#'
#' @examples
#'   NULL
adjust_and_sort <- function(pVals_vec,
                            genesets_ls,
                            adjust = TRUE,
                            proc_vec = c("Bonferroni",
                                         "Holm",
                                         "Hochberg",
                                         "SidakSS",
                                         "SidakSD",
                                         "BH",
                                         "BY",
                                         "ABH",
                                         "TSBH"),
                            ...){

  pVals_df <- data.frame(pathways = names(pVals_vec),
                         setsize = genesets_ls$setsize,
                         terms = genesets_ls$TERMS,
                         rawp = unname(pVals_vec),
                         stringsAsFactors = FALSE)
  rownames(pVals_df) <- NULL

  if(adjust){

    adjusted_mat <- adjustRaw_pVals(rawp = pVals_df$rawp,
                                    proc = proc_vec,
                                    ...)
    adjusted_mat <- adjusted_mat[, -1, drop = FALSE]

    pVals_df <- cbind(pVals_df, adjusted_mat)
    out_df <- pVals_df[order(pVals_df[, proc_vec[1]], pVals_df$rawp), ]

  } else {
    out_df <- pVals_df[order(pVals_df$rawp), ]
  }

  out_df

}
