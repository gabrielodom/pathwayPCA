#' AES-PCA permutation test of continuous response for pathway PCs
#'
#' @description Given an \code{OmicsReg} object and a list of pathway PCs from
#'   the \code{\link{extract_aesPCs}} function, test if each expressed pathway
#'   in the bio-assay design matrix is significantly related to the continuous
#'   response.
#'
#' @param OmicsReg A data object of class \code{OmicsReg}, created by the
#'   \code{\link{create_OmicsReg}} function.
#' @param pathwayPCs_ls A list of pathway PC matrices returned by the
#'   \code{\link{extract_aesPCs}} function.
#' @param numReps How many permuted models to fit? Defaults to 1000.
#' @param parallel Should the computation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation? Defaults to \code{NULL}.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return A named vector of pathway permutation \eqn{p}-values.
#'
#' @details This function takes in a list of the first principal components
#'   from each pathway and an object of class \code{OmicsReg}. This function
#'   will then calculate the AIC of a multivariate linear model (via the
#'   \code{\link[stats]{lm}} function) with the original observations as
#'   response and the pathway principal components as the predictor matrix.
#'
#'   Then, this function will create \code{numReps} permutations of the
#'   regression response, fit models to each of these premuted responses
#'   (holding the path predictor matrix fixed), and calculate the AIC of each
#'   model. This function will return a named vector of permutation
#'   \eqn{p}-values, where the value for each pathway is the proportion of
#'   models for which the AIC of the permuted response model is less than the
#'   AIC of the original model.
#'
#' @seealso \code{\link{create_OmicsReg}}; \code{\link{extract_aesPCs}};
#'   \code{\link[stats]{lm}}; \code{\link{sample_Regresp}}
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsReg.R
#' @include aesPC_extract_OmicsPath_PCs.R
#'
#' @importFrom methods setGeneric
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
#'
#' @rdname permTest_OmicsReg
setGeneric("permTest_OmicsReg",
           function(OmicsReg,
                    pathwayPCs_ls,
                    numReps = 1000,
                    parallel = FALSE,
                    numCores = NULL,
                    ...){
             standardGeneric("permTest_OmicsReg")
           }
)

#' @importFrom stats lm
#' @importFrom stats AIC
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#'
#' @rdname permTest_OmicsReg
setMethod(f = "permTest_OmicsReg", signature = "OmicsReg",
          definition = function(OmicsReg,
                                pathwayPCs_ls,
                                numReps = 1000,
                                parallel = FALSE,
                                numCores = NULL,
                                ...){

            ###  Function Setup  ###
            permute_RegFit <- function(pathwayPCs_mat,
                                       obj_OmicsReg,
                                       numReps_int = numReps,
                                       parametric = FALSE){
              # browser()

              ###  True Model  ###
              response <- obj_OmicsReg@response
              trueAIC <- AIC(lm(response ~ pathwayPCs_mat))


              ###  Permuted Model  ###
              permuteAIC_fun <- function(){

                perm_resp <- sample_Regresp(obj_OmicsReg@response,
                                            parametric = parametric)
                AIC(lm(perm_resp ~ pathwayPCs_mat))

              }

              permAIC <- replicate(n = numReps_int, expr = permuteAIC_fun())

              ###  Return  ###
              mean(permAIC < trueAIC)

            }

            ###  Computation  ###
            # browser()

            if(parallel){
              # browser()

              ###  Parallel Computing Setup  ###
              message("Initializing Cluster")
              # require(parallel)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(deparse(quote(OmicsReg)),
                                 deparse(quote(numReps)))
              clusterExport(cl = clust,
                            varlist = clustVars_vec,
                            envir = environment())
              invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
              message("DONE")

              ###  Extract PCs  ###
              message("Extracting Pathway p-Values in Parallel")
              pValues_vec <- parSapply(cl = clust,
                                       pathwayPCs_ls,
                                       permute_RegFit,
                                       obj_OmicsReg = OmicsReg,
                                       numReps_int = numReps)
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway p-Values Serially")
              pValues_vec <- sapply(pathwayPCs_ls,
                                    permute_RegFit,
                                    obj_OmicsReg = OmicsReg,
                                    numReps_int = numReps)
              message("DONE")

            }

            pValues_vec

          })
