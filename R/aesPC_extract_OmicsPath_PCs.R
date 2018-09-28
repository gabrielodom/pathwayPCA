#' Extract AES-PCs from recorded pathway-subsets of a mass spectrometry or
#'   bio-assay data frame
#'
#' @description Given a clean \code{OmicsPath} object (cleaned by the
#'   \code{\link{IntersectOmicsPwyCollct}} function), extract the first
#'   principal components (PCs) from each pathway with features recorded in the
#'   assay design matrix.
#'
#' @param object An object of class \code{OmicsPathway}.
#' @param numPCs The number of PCs to extract from each pathway. Defaults to 1.
#' @param parallel Should the computation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation? Defaults to \code{NULL}.
#' @param standardPCA Should the function return the AES-PCA PCs and loadings
#'   (\code{FALSE}) or the standard PCA PCs and loadings (\code{TRUE})? Defaults
#'   to \code{FALSE}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return Two lists of matrices: \code{PCs} and \code{loadings}. Each element
#'   of both lists will be named by its pathway. The elements of the \code{PCs}
#'   list will be \eqn{N \times} \code{numPCs} matrices containing the first
#'   \code{numPCs} principal components from each pathway. The elements of the
#'   \code{loadings} list will be \code{numPCs} \eqn{\times p} projection
#'   matrices containing the loadings corresponding to the first \code{numPCs}
#'   principal components from each pathway. See "Details" for more information.
#'
#' @details This function takes in a data frame with named columns and a pathway
#'   list as an \code{OmicsPathway} object which has had unrecorded -Omes
#'   removed from the corresponding pathway collection by the
#'   \code{\link{IntersectOmicsPwyCollct}} function. This function will then
#'   iterate over the list of pathways, extracting columns from the assay design
#'   matrix which match the genes listed in that pathway as a sub-matrix (as a
#'   \code{data.frame} object). This function will then call the
#'   \code{\link{aespca}} on each data frame in the list of pathway-specific
#'   design matrices, extracting the first \code{numPCs} AES principal
#'   components from each pathway data frame. These PC matrices are returned as
#'   a named list.
#'
#'   NOTE: some genes will be included in more than one pathway, so these
#'   pathways are not mutually exclusive. Further note that there may be many
#'   genes in the assay design matrix that are not included in the pathways,
#'   so these will not be extracted to the list. It is then vitally important to
#'   use either a very broad and generic list of pathways or a pathways list
#'   that is compatible to the assay data supplied.
#'
#' @seealso \code{\link{CreateOmicsPath}}; \code{\link{aespca}};
#'    \code{\link{IntersectOmicsPwyCollct}}
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @importFrom methods setGeneric
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
#'
#' @rdname ExtractAESPCs
setGeneric("ExtractAESPCs",
           function(object, numPCs = 1,
                    parallel = FALSE, numCores = NULL,
                    standardPCA = FALSE,
                    ...){
             standardGeneric("ExtractAESPCs")
           }
)

#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parLapplyLB
#' @importFrom parallel stopCluster
#'
#' @rdname ExtractAESPCs
setMethod(f = "ExtractAESPCs", signature = "OmicsPathway",
          definition = function(object,
                                numPCs = 1,
                                parallel = FALSE,
                                numCores = NULL,
                                standardPCA = FALSE,
                                ...){
            # browser()
            pathSets_ls <- object@trimPathwayCollection
            data_Omes <- lapply(pathSets_ls$pathways, function(x){
              object@assayData_df[x]
            })

            if(parallel){
              # browser()

              ###  Parallel Computing Setup  ###
              message("Initializing Computing Cluster: ", appendLF = FALSE)
              # require(parallel)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(deparse(quote(data_Omes)),
                                 deparse(quote(numPCs)))
              clusterExport(cl = clust,
                            varlist = clustVars_vec,
                            envir = environment())
              invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
              message("DONE")

              ###  Extract PCs  ###
              message("Extracting Pathway PCs in Parallel: ", appendLF = FALSE)
              out_ls <- parLapplyLB(cl = clust,
                                    data_Omes,
                                    function(pathway_df){
                                      aespca(X = pathway_df,
                                             d = numPCs)
                                    })
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway PCs Serially: ", appendLF = FALSE)
              out_ls <- lapply(data_Omes,
                               function(path_df){
                                 aespca(X = path_df,
                                        d = numPCs)
                               })
              message("DONE")

            }


            ###  Return  ###
            if(standardPCA){

              PCs_ls      <- lapply(out_ls, `[[`, "oldScore")
              loadings_ls <- lapply(out_ls, `[[`, "oldLoad")

            } else {

              PCs_ls      <- lapply(out_ls, `[[`, "aesScore")
              loadings_ls <- lapply(out_ls, `[[`, "aesLoad")

            }

            list(
              PCs = PCs_ls,
              loadings = loadings_ls
            )

          })
