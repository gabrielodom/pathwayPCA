#' AES-PCA permutation test of survival response for pathway PCs
#'
#' @description Given an \code{OmicsSurv} object and a list of pathway principal
#'    components (PCs) from the \code{\link{ExtractAESPCs}} function, test if
#'    each pathway with features recorded in the bio-assay design matrix is
#'    significantly related to the survival output.
#'
#' @param OmicsSurv A data object of class \code{OmicsSurv}, created by the
#'   \code{\link{CreateOmicsSurv}} function.
#' @param pathwayPCs_ls A list of pathway PC matrices returned by the
#'   \code{\link{ExtractAESPCs}} function.
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
#'   from each pathway and an object of class \code{OmicsSurv}. This function
#'   will then calculate the AIC of a Cox Proportional Hazards model (via the
#'   \code{\link[survival]{coxph}} function) with the original observations as
#'   response and the pathway principal components as the predictor matrix.
#'
#'   Then, this function will create \code{numReps} permutations of the survival
#'   response, fit models to each of these permuted responses (holding the path
#'   predictor matrix fixed), and calculate the AIC of each model. This function
#'   will return a named vector of permutation \eqn{p}-values, where the value
#'   for each pathway is the proportion of models for which the AIC of the
#'   permuted response model is less than the AIC of the original model.
#'
#' @seealso \code{\link{CreateOmicsSurv}}; \code{\link{ExtractAESPCs}};
#'   \code{\link[survival]{coxph}}; \code{\link{sample_Survivalresp}}
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsSurv.R
#' @include aesPC_extract_OmicsPath_PCs.R
#'
#' @importFrom methods setGeneric
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
#'
#' @rdname permTest_OmicsSurv
setGeneric("permTest_OmicsSurv",
           function(OmicsSurv,
                    pathwayPCs_ls,
                    numReps = 1000,
                    parallel = FALSE,
                    numCores = NULL,
                    ...){
             standardGeneric("permTest_OmicsSurv")
           }
)

#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom stats AIC
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#'
#' @rdname permTest_OmicsSurv
setMethod(f = "permTest_OmicsSurv", signature = "OmicsSurv",
          definition = function(OmicsSurv,
                                pathwayPCs_ls,
                                numReps = 1000,
                                parallel = FALSE,
                                numCores = NULL,
                                ...){

            ###  Function Setup  ###
            permute_SurvFit <- function(pathwayPCs_mat,
                                        obj_OmicsSurv,
                                        numReps_int = numReps,
                                        parametric = FALSE){
              # browser()

              ###  True Model  ###
              # I know I'm doing this more than once, but I don't know how to
              #   rewrite the sample_Survivalresp() function to take in a
              #   Surv object
              response <- Surv(time = obj_OmicsSurv@eventTime,
                               event = obj_OmicsSurv@eventObserved)
              trueAIC <- AIC(coxph(response ~ pathwayPCs_mat))


              ###  Permuted Model  ###
              permuteAIC_fun <- function(){

                perm_resp <- sample_Survivalresp(response_vec = obj_OmicsSurv@eventTime,
                                                 event_vec = obj_OmicsSurv@eventObserved,
                                                 parametric = parametric)
                perm_Surv <- Surv(time = perm_resp$response_vec,
                                  event = perm_resp$event_vec)
                AIC(coxph(perm_Surv ~ pathwayPCs_mat))

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
              message("Initializing Computing Cluster: ", appendLF = FALSE)
              # require(parallel)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(deparse(quote(OmicsSurv)),
                                 deparse(quote(numReps)))
              clusterExport(cl = clust,
                            varlist = clustVars_vec,
                            envir = environment())
              invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
              message("DONE")

              ###  Extract PCs  ###
              message("Extracting Pathway p-Values in Parallel: ",
                      appendLF = FALSE)
              pValues_vec <- parSapply(cl = clust,
                                       pathwayPCs_ls,
                                       permute_SurvFit,
                                       obj_OmicsSurv = OmicsSurv,
                                       numReps_int = numReps)
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway p-Values Serially")
              pValues_vec <- sapply(pathwayPCs_ls,
                                    permute_SurvFit,
                                    obj_OmicsSurv = OmicsSurv,
                                    numReps_int = numReps)
              message("DONE")

            }

            pValues_vec

          })
