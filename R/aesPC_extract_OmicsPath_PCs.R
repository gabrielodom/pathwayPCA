#' Extract AES-PCs from expressed pathway-subsets of a mass spectrometry or
#'   bio-assay data frame
#'
#' @description Given a clean \code{OmicsPath} object (cleaned by the
#'   \code{\link{expressedOmes}} function), extract the first principal
#'   components from each expressed pathway in the assay design matrix.
#'
#' @param object An object of class \code{OmicsPathway}.
#' @param trim The minimum cutoff of expressed -Ome measures before a pathway
#'   is excluded. Defaults to 3.
#' @param numPCs The number of PCs to extract from each pathway. Defaults to 1.
#' @param parallel Should the computation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation? Defaults to \code{NULL}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return A list of matrices. Each element of the list will be named
#'   by its pathway, and the elements will be \eqn{N \times} \code{numPCs}
#'   matrices containing the first \code{numPCs} principal components from each
#'   pathway. See "Details" for more information.
#'
#' @details This function takes in a data frame with named columns and a pathway
#'   list as an \code{OmicsPathway} object which has had unexpressed -Omes
#'   removed by the \code{\link{expressedOmes}} function. This function will
#'   then iterate over the list of pathways, extracting columns from the assay
#'   design matrix which match the genes listed in that pathway as a sub-matrix
#'   (as a \code{data.frame} object). This function will then call the
#'   \code{\link{aespca}} on each data frame in the list of pathway-specific
#'   design matrices, extracting the first \code{numPCs} AES principal
#'   components from each pathway data frame. These PC matrices are returned as
#'   a named list.
#'
#'   NOTE: some genes will be included in more than one pathway, so these
#'   pathways are not mutually exclusive. Further note that there may be many
#'   genes in the assay design matrix that are not included in the pathway sets,
#'   so these will not be extracted to the list. It is then vitally important to
#'   use either a very broad and generic pathway set list or a pathway set list
#'   that is appropriate for the assay data supplied.
#'
#' @seealso \code{\link{create_OmicsPath}}; \code{\link{expressedOmes}};
#'   \code{\link{aespca}}
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#' @include subsetExpressed-omes.R
#'
#' @importFrom methods setGeneric
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
#'
#' @rdname extract_aesPCs
setGeneric("extract_aesPCs",
           function(object, trim = 3, numPCs = 1,
                    parallel = FALSE, numCores = NULL,
                    ...){
             standardGeneric("extract_aesPCs")
           }
)

#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parLapplyLB
#' @importFrom parallel stopCluster
#'
#' @rdname extract_aesPCs
setMethod(f = "extract_aesPCs", signature = "OmicsPathway",
          definition = function(object,
                                trim = 3,
                                numPCs = 1,
                                parallel = FALSE,
                                numCores = NULL,
                                ...){
            # browser()
            pathSets_ls <- object@pathwaySet
            data_Omes <- lapply(pathSets_ls$pathways, function(x){
              object@assayData_df[x]
            })

            if(parallel){
              # browser()

              ###  Parallel Computing Setup  ###
              message("Initializing Computing Cluster")
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
              message("Extracting Pathway PCs in Parallel")
              PCs_ls <- parLapplyLB(cl = clust,
                                    data_Omes,
                                    function(pathway_df){
                                      aespca(X = pathway_df,
                                             d = numPCs)$score
                                    })
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway PCs Serially")
              PCs_ls <- lapply(data_Omes,
                               function(path_df){
                                 aespca(X = path_df,
                                        d = numPCs)$score
                               })
              message("DONE")

            }

            PCs_ls

          })
