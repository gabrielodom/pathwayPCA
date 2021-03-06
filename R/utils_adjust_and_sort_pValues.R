#' Tabulate, adjust, and sort pathway \eqn{p}-values
#'
#' @description Adjust the pathway \eqn{p}-values, then return a data frame of
#'   the relevant pathway information, sorted by adjusted significance.
#'
#' @param pVals_vec A named vector of permutation \eqn{p}-values returned by the
#'   \code{\link{PermTestSurv}}, \code{\link{PermTestReg}}, or
#'   \code{\link{PermTestCateg}} functions when the analysis performed was
#'   AES-PCA. Otherwise, when the analysis was performed with Supervised PCA, a
#'   named vector of \eqn{p}-values from the \code{\link{GumbelMixpValues}}
#'   function.
#' @param genesets_ls A list of known gene pathways, trimmed to match the given
#'   assay data by the \code{\link{IntersectOmicsPwyCollct}} function. This
#'   pathway list must contain:
#'   \itemize{
#'     \item{\code{pathways} : }{A named list of character vectors where each
#'        vector contains the names of the genes in that specific pathway.}
#'     \item{\code{TERMS} : }{A character vector the same length as
#'        \code{pathways} containing the full pathway descriptions.}
#'     \item{\code{n_tested} : }{An integer vector the same length as
#'        \code{pathways} containing the number of genes present in the pathway
#'        after trimming. Pathways list trimming is done in the
#'       \code{\link{IntersectOmicsPwyCollct}} function.}
#'   }
#' @param adjust Should you adjust the \eqn{p}-values for multiple comparisons?
#'   Defaults to TRUE.
#' @param proc_vec Character vector of procedures. The returned data frame will
#'   be sorted in ascending order by the first procedure in this vector, with
#'   ties broken by the unadjusted \eqn{p}-value. If only one procedure is
#'   selected, then it is necessarily the first procedure. Defaults to
#'   \code{"BH"} (Benjamini and Hochberg, 1995).
#' @param ... Additional arguments to pass to the \code{\link{ControlFDR}}
#'   function.
#'
#' @return A data frame with columns
#' \itemize{
#'   \item{\code{pathways} : }{The names of the pathways in the \code{Omics*}}
#'     object (stored in \code{object@@trimPathwayCollection$pathways}).
#'   \item{\code{n_tested} : }{The number of genes in each pathway after being
#'    trimmed to match the assay. Given in the \code{n_tested} element of the
#'    trimmed pathway collection.}
#'   \item{\code{terms} : }{The pathway title, as stored in the
#'     \code{object@@trimPathwayCollection$TERMS} object.}
#'   \item{\code{description} : }{The pathway description, if it is stored in
#'     the \code{object@@trimPathwayCollection$description} object.}
#'   \item{\code{rawp} : }{The unadjusted \eqn{p}-values of each pathway.}
#'   \item{\code{...} : }{Additional columns as specified through the
#'     \code{adjustment} argument.}
#' }
#'
#' The data frame will be sorted in ascending order by the method specified
#'   first in the \code{adjustment} argument. If \code{adjustpValues = FALSE},
#'   then the data frame will be sorted by the raw \eqn{p}-values. If you have
#'   the suggested \code{tidyverse} package suite loaded, then this data frame
#'   will print as a \code{\link[tibble]{tibble}}. Otherwise, it will stay a
#'   simple data frame.
#'
#' @details This is a wrapper function for the \code{\link{ControlFDR}}
#'   function. The number of \eqn{p}-values passed to the \code{pVals_vec}
#'   argument \emph{must} equal the number of pathways and set size values in
#'   the \code{genesets_ls} argument. If you trimmed a pathway from \eqn{p}-
#'   value calculation, then pad this missing value with an \code{NA}.
#'
#' @keywords internal
#'
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Call this function through AESPCA_pVals() or SuperPCA_pVals() instead.
#'
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_Omics <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(2:3)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:3],
#'     respType = "surv"
#'   )
#'
#'   ###  Extract Pathway PCs and Loadings  ###
#'   colonPCs_ls <- ExtractAESPCs(
#'     object = colon_Omics,
#'     parallel = TRUE,
#'     numCores = 2
#'   )
#'
#'   ###  Pathway p-Values  ###
#'   pVals <- PermTestSurv(
#'     OmicsSurv = colon_Omics,
#'     pathwayPCs_ls = colonPCs_ls$PCs,
#'     parallel = TRUE,
#'     numCores = 2
#'   )
#'
#'   ###  Create Table of p-Values  ###
#'   trimmed_PC <- getTrimPathwayCollection(colon_Omics)
#'   TabulatepValues(
#'     pVals_vec = pVals,
#'     genesets_ls = trimmed_PC
#'   )
#' }
#'
TabulatepValues <- function(pVals_vec,
                            genesets_ls,
                            adjust = TRUE,
                            proc_vec = c("BH",
                                         "Bonferroni",
                                         "Holm",
                                         "Hochberg",
                                         "SidakSS",
                                         "SidakSD",
                                         "BY",
                                         "ABH",
                                         "TSBH"),
                            ...){

  # browser()

  # # If we remove a pathway, then the pValue for that pathway should
  # #   exist, and it should be NA. Right now, if we remove 10 pathways,
  # #   then the pVals_vec will be 10 shorter. We need these vectors to
  # #   be the same size.
  # pValsFull_vec <- rep(NA, length(genesets_ls$TERMS))
  # names(pValsFull_vec) <- names(genesets_ls$TERMS)
  # pValsFull_vec[names(pVals_vec)] <- pVals_vec
  # # We are passing in the trimmed pathway collection anyway, so this is now
  # #   pointless.

  if(is.null(genesets_ls$description)){
    pVals_df <- data.frame(
      pathways = names(genesets_ls$TERMS),
      n_tested = genesets_ls$n_tested,
      terms = genesets_ls$TERMS,
      rawp = unname(pVals_vec),
      stringsAsFactors = FALSE
    )
  } else {
    pVals_df <- data.frame(
      pathways = names(genesets_ls$TERMS),
      n_tested = genesets_ls$n_tested,
      terms = genesets_ls$TERMS,
      description = genesets_ls$description,
      rawp = unname(pVals_vec),
      stringsAsFactors = FALSE
    )
  }
  rownames(pVals_df) <- NULL

  if(adjust){

    proc_vec <- match.arg(proc_vec, several.ok = TRUE)
    proc_vec <- unique(proc_vec)

    adjusted_df <- ControlFDR(
      rawp = pVals_df$rawp,
      proc = proc_vec,
      na.rm = anyNA(pVals_df$rawp),
      ...
    )

    TSBHplace_idx <- which(proc_vec == "TSBH")

    if(length(TSBHplace_idx) == 0){

      # If the "TSBH" method isn't included
      orderedNames <- proc_vec

    } else {

      # if the "TSBH" procedure is included
      TSBHnames_idx <- grep("TSBH", colnames(adjusted_df))
      TSBHnames <- colnames(adjusted_df)[TSBHnames_idx]
      nProcs <- length(proc_vec)

      if(TSBHplace_idx == 1){

        # "TSBH" is the first procedure
        orderedNames <- c(TSBHnames, proc_vec[-1])

      } else if(TSBHplace_idx == nProcs) {

        # "TSBH" is the last procedure
        orderedNames <- c(proc_vec[-nProcs], TSBHnames)

      } else {

        # "TSBH" is some procedure in the middle of the proc list
        orderedNames <- c(
          proc_vec[seq_len(TSBHplace_idx - 1)],
          TSBHnames,
          proc_vec[(TSBHplace_idx + 1):nProcs]
        )

      }

    }

    # browser()

    adjusted_df <- adjusted_df[, orderedNames, drop = FALSE]

    fwerNames <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD")
    newAdjNames <- vapply(
      orderedNames,
      function(i){
        ifelse(
          i %in% fwerNames,
          paste0("FWER_", i),
          paste0("FDR_", i)
        )
      },
      FUN.VALUE = character(1)
    )
    colnames(adjusted_df) <- newAdjNames

    pVals_df <- cbind(pVals_df, adjusted_df)
    out_df <- pVals_df[order(pVals_df[, newAdjNames[1]], pVals_df$rawp), ]

  } else {
    out_df <- pVals_df[order(pVals_df$rawp), ]
  }

  class(out_df) <- c("tbl_df", "tbl", "data.frame")
  out_df

}
