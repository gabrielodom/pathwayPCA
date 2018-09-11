#' Test pathway association with AES-PCA
#'
#' @description Given a supervised \code{OmicsPath} object (one of
#'    \code{OmicsSurv}, \code{OmicsReg}, or \code{OmicsCateg}), extract the
#'    first \eqn{k} adaptive, elastic-net, sparse principal components (PCs)
#'    from each expressed pathway in the -Omics assay design matrix, test their
#'    association with the response matrix, and return a data frame of the
#'    adjusted \eqn{p}-values for each pathway.
#'
#' @param object An object of class \code{OmicsPathway} with a response matrix
#'   or vector.
#' @param numPCs The number of PCs to extract from each pathway. Defaults to 1.
#' @param numReps The number of permutations to take of the data to calculate a
#'   \eqn{p}-value for each pathway. Defaults to 1000.
#' @param parallel Should the computation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation? Defaults to \code{NULL}.
#' @param asPCA Should the computation return the eigenvectors and eigenvalues
#'   instead of the adaptive, elastic-net, sparse principal components and their
#'   corresponding loadings. Defaults to \code{FALSE}; this should be used for
#'   diagnostic or comparative purposes only.
#' @param adjustpValues Should you adjust the \eqn{p}-values for multiple
#'   comparisons? Defaults to TRUE.
#' @param adjustment Character vector of procedures. The returned data frame
#'   will be sorted in ascending order by the first procedure in this vector,
#'   with ties broken by the unadjusted \eqn{p}-value. If only one procedure is
#'   selected, then it is necessarily the first procedure. See the documentation
#'   for the \code{\link{adjustRaw_pVals}} function for the adjustment procedure
#'   definitions and citations.
#' @param ... Dots for additional internal arguments.
#'
#' @return A results list with class \code{aespcOut}. This list has three
#'    components: a data frame of pathway details, pathway \eqn{p}-values, and
#'    potential adjustments to those values (\code{pVals_df}); a list of the
#'    first \code{numPCs} \emph{score} vectors for each pathway (\code{PCs_ls});
#'    and a list of the first \code{numPCs} feature loading vectors for each
#'    pathway (\code{loadings_ls}). The \eqn{p}-value data frame has columns:
#' \itemize{
#'   \item{\code{pathways} : }{The names of the pathways in the \code{Omics*}}
#'     object (given in \code{object@@trimPathwayCollection$pathways}.)
#'   \item{\code{setsize} : }{The number of genes in each of the original
#'     pathways (given in the \code{object@@trimPathwayCollection$setsize}
#'     object).}
#'   \item{\code{trim_size} : }{The number of genes in each of the trimmed
#'     pathways (given in the \code{object@@trimPathwayCollection$trim_setsize}
#'     object).}
#'   \item{\code{terms} : }{The pathway description, as given in the
#'     \code{object@@trimPathwayCollection$TERMS} object.}
#'   \item{\code{rawp} : }{The unadjusted \eqn{p}-values of each pathway.}
#'   \item{\code{...} : }{Additional columns of adjusted \eqn{p}-values as
#'     specified through the \code{adjustment} argument.}
#' }
#'
#' The data frame will be sorted in ascending order by the method specified
#'    first in the \code{adjustment} argument. If \code{adjustpValues = FALSE},
#'    then the data frame will be sorted by the raw \eqn{p}-values. If you have
#'    the suggested \code{tidyverse} package suite loaded, then this data frame
#'    will print as a \code{\link[tibble]{tibble}}. Otherwise, it will print as
#'    a data frame.
#'
#' @details This is a wrapper function for the \code{\link{extract_aesPCs}},
#'    \code{\link{permTest_OmicsSurv}}, \code{\link{permTest_OmicsReg}}, and
#'    \code{\link{permTest_OmicsCateg}} functions.
#'
#'   Please see our Quickstart Guide for this package:
#'   \url{https://gabrielodom.github.io/pathwayPCA/articles/C1-Quickstart_Guide.html}
#'
#' @seealso \code{\link{create_OmicsPath}}; \code{\link{create_OmicsSurv}};
#'    \code{\link{create_OmicsReg}}; \code{\link{create_OmicsCateg}};
#'    \code{\link{extract_aesPCs}}; \code{\link{permTest_OmicsSurv}};
#'    \code{\link{permTest_OmicsReg}}; \code{\link{permTest_OmicsCateg}};
#'    \code{\link{adjust_and_sort}}
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#' @include createClass_OmicsSurv.R
#' @include createClass_OmicsReg.R
#' @include createClass_OmicsCateg.R
#' @include aesPC_permtest_CoxPH.R
#' @include aesPC_permtest_LS.R
#' @include aesPC_permtest_GLM.R
#'
#' @importFrom methods setGeneric
#'
#' @examples
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_OmicsSurv(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     eventTime_num = colonSurv_df$OS_time,
#'     eventObserved_lgl = as.logical(colonSurv_df$OS_event)
#'   )
#'
#'   ###  Calculate Pathway p-Values  ###
#'   colonSurv_pVals_df <- AESPCA_pVals(
#'     object = colon_OmicsSurv,
#'     numReps = 5000,
#'     parallel = TRUE,
#'     numCores = 16,
#'     adjustpValues = TRUE,
#'     adjustment = c("Hoch", "SidakSD")
#'   )
#' }
#'
#' @rdname AESPCA_pVals
setGeneric("AESPCA_pVals",
           function(object,
                    numPCs = 1,
                    numReps = 1000,
                    parallel = FALSE,
                    numCores = NULL,
                    asPCA = FALSE,
                    adjustpValues = TRUE,
                    adjustment = c("Bonferroni",
                                   "Holm",
                                   "Hochberg",
                                   "SidakSS",
                                   "SidakSD",
                                   "BH",
                                   "BY",
                                   "ABH",
                                   "TSBH"),
                    ...){
             standardGeneric("AESPCA_pVals")
           }
)

#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#'
#' @rdname AESPCA_pVals
setMethod(f = "AESPCA_pVals", signature = "OmicsPathway",
          definition = function(object,
                                numPCs = 1,
                                numReps = 1000,
                                parallel = FALSE,
                                numCores = NULL,
                                asPCA = FALSE,
                                adjustpValues = TRUE,
                                adjustment = c("Bonferroni",
                                               "Holm",
                                               "Hochberg",
                                               "SidakSS",
                                               "SidakSD",
                                               "BH",
                                               "BY",
                                               "ABH",
                                               "TSBH"),
                                ...){
            # browser()

            ###  Calculate AES-PCs  ###
            message(" Part 1: Calculate Pathway AES-PCs\n")
            aespca_ls <- extract_aesPCs(
              object = object,
              numPCs = numPCs,
              parallel = parallel,
              numCores = numCores
            )


            if(asPCA){

              PCs_ls      <- lapply(aespca_ls, `[[`, "oldScore")
              loadings_ls <- lapply(aespca_ls, `[[`, "oldLoad")

            } else {

              PCs_ls      <- lapply(aespca_ls, `[[`, "aesScore")
              loadings_ls <- lapply(aespca_ls, `[[`, "aesLoad")

            }
            rm(aespca_ls)


            ###  Permutation Pathway p-Values  ###
            message("\n Part 2: Calculate Permuted Pathway p-Values\n")
            obj_class <- class(object)
            switch(obj_class,
                   OmicsSurv = {
                     pVals_vec <- permTest_OmicsSurv(OmicsSurv = object,
                                                     pathwayPCs_ls = PCs_ls,
                                                     numReps = numReps,
                                                     parallel = parallel,
                                                     numCores = numCores)
                   },
                   OmicsReg = {
                     pVals_vec <- permTest_OmicsReg(OmicsReg = object,
                                                    pathwayPCs_ls = PCs_ls,
                                                    numReps = numReps,
                                                    parallel = parallel,
                                                    numCores = numCores)
                   },
                   OmicsCateg = {
                     pVals_vec <- permTest_OmicsCateg(OmicsCateg = object,
                                                      pathwayPCs_ls = PCs_ls,
                                                      numReps = numReps,
                                                      parallel = parallel,
                                                      numCores = numCores)
                   }
            )


            ###  Adjust Pathway p-Values  ###
            if(adjustpValues){

              message("\n Part 3: Adjusting p-Values and Sorting Pathway p-Value Data Frame\n")
              adjustment <- match.arg(adjustment, several.ok = TRUE)

            } else {
              message("\n Part 3: Sorting Pathway p-Value Data Frame\n")
            }

            pathwayGeneSets_ls <- object@trimPathwayCollection
            out_df <- adjust_and_sort(pVals_vec = pVals_vec,
                                      genesets_ls = pathwayGeneSets_ls,
                                      adjust = adjustpValues,
                                      proc_vec = adjustment,
                                      ...)
            message("DONE")


            ###  Re-order PCA Output  ###
            PCs_ls <- PCs_ls[out_df$pathways]
            loadings_ls <- loadings_ls[out_df$pathways]


            ###  Return  ###
            out_ls <- list(
              pVals_df    = out_df,
              PCs_ls      = PCs_ls,
              loadings_ls = loadings_ls
            )

            class(out_ls) <- c("aespcOut", "pathwayPCA", "list")
            out_ls

          })
