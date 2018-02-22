#' Extract PCs from pathway-subsets of a mass spectrometry matrix
#'
#' @description Given an \code{OmicsPath} object, extract the first principal
#'   components from each expressed pathway in the MS design matrix
#'
#' @param object An object of class \code{OmicsPathway}
#' @param trim The minimum cutoff of expressed -ome measures before a pathway
#'   is excluded. Defaults to 3.
#' @param numPCs The number of PCs to extract from each pathway. Defaults to 1.
#' @param parallel Should the comuptation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation?
#' @param ... Dots for additional internal arguments (currently unused)
#'
#' @return A list of pathway PC matrices. Each element of the list will be named
#'   by its pathway, and the elements will be \eqn{N \times} \code{numPCs}
#'   matrices containing the first \code{numPCs} principal components from each
#'   pathway. See "details" for more information
#'
#' @details This function takes in a data frame with named columns and a pathway
#'   list as an \code{OmicsPathway} object. This function will then call the
#'   \code{\link{expressedOmes}} function to iterate over the list of pathways,
#'   extracting columns from the MS design matrix which match the genes listed
#'   in that pathway and removing any pathways with fewer than \code{trim}
#'   expressed genes. This function will then call the \code{\link{aespca}}
#'   on the list of pathway-specific data frames returned by the
#'   \code{\link{expressedOmes}} function, extracting the first \code{numPCs}
#'   from each pathway data frame. These PC matrices are returned as a named
#'   list.
#'
#'   NOTE: some genes will be included in more than one pathway, so these
#'   pathways are not a partition of the MS matrix supplied. Further, note that
#'   there may be many genes in the MS design matrix that are not included in
#'   the pathway sets, so these will not be extracted to the list. It is then
#'   vitally important to use either a very broad and generic pathway set list
#'   or a pathway set list that is appropriate for the mass spectrometry data
#'   supplied.
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

            data_Omes <- expressedOmes(object, returnClass = "list", trim = trim)
            # data_Omes <- data_Omes[1:18]
            n <- nrow(data_Omes[[1]])

            if(parallel){
              # browser()

              ###  Parallel Computing Setup  ###
              message("Initializing Cluster")
              # require(parallel)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(deparse(quote(data_Omes)),
                                 deparse(quote(n)),
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
                                             n = n, d = numPCs,
                                             type = "predictor")$score
                                    })
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway PCs Serially")
              PCs_ls <- lapply(data_Omes,
                               function(path_df){
                                 aespca(X = path_df,
                                        n = n, d = numPCs,
                                        type = "predictor")$score
                               })
              message("DONE")

            }

            PCs_ls

          })
