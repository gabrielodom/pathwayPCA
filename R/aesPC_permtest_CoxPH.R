#' AES-PCA permutation test of survival response for pathway PCs
#'
#' @description Given an \code{OmicsSurv} object and a list of pathway principal
#'    components (PCs) from the \code{\link{ExtractAESPCs}} function, test if
#'    each pathway with features recorded in the bio-assay design matrix is
#'    significantly related to the survival output.
#'
#' @param OmicsSurv A data object of class \code{OmicsSurv}, created by the
#'   \code{\link{CreateOmics}} function.
#' @param pathwayPCs_ls A list of pathway PC matrices returned by the
#'   \code{\link{ExtractAESPCs}} function.
#' @param numReps How many permutations to estimate the \eqn{p}-value? Defaults
#'    to 0 (that is, to estimate the \eqn{p}-value parametrically). If
#'    \code{numReps} > 0, then the non-parametric, permutation \eqn{p}-value
#'    will be returned based on the number of random samples specified.
#' @param parallel Should the computation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation? Internally defaults to the number of available cores minus 2.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return A named vector of pathway permutation \eqn{p}-values.
#'
#' @details This function takes in a list of the first principal components
#'   from each pathway and an object of class \code{OmicsSurv}. This function
#'   will then calculate the AIC of a Cox Proportional Hazards model (via the
#'   \code{\link[survival]{coxph}} function) with the original observations as
#'   response and the pathway principal components as the predictor matrix. Note
#'   that the AIC and log-likelihood are proportional because the number of
#'   parameters in each pathway is constant.
#'
#'   Then, this function will create \code{numReps} permutations of the survival
#'   response, fit models to each of these permuted responses (holding the path
#'   predictor matrix fixed), and calculate the AIC of each model. This function
#'   will return a named vector of permutation \eqn{p}-values, where the value
#'   for each pathway is the proportion of models for which the AIC of the
#'   permuted response model is less than the AIC of the original model.
#'
#' @seealso \code{\link{CreateOmics}}; \code{\link{ExtractAESPCs}};
#'   \code{\link[survival]{coxph}}; \code{\link{SampleSurv}}
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsSurv.R
#' @include aesPC_extract_OmicsPath_PCs.R
#'
#' @importFrom methods setGeneric
#'
#' @keywords internal
#'
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTION DIRECTLY.
#'   # Use AESPCA_pVals() instead
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
#'   PermTestSurv(
#'     OmicsSurv = colon_Omics,
#'     pathwayPCs_ls = colonPCs_ls$PCs,
#'     parallel = TRUE,
#'     numCores = 2
#'   )
#' }
#'
#' @rdname PermTestSurv
setGeneric("PermTestSurv",
           function(OmicsSurv,
                    pathwayPCs_ls,
                    numReps = 0L,
                    parallel = FALSE,
                    numCores = NULL,
                    ...){
             standardGeneric("PermTestSurv")
           }
)

#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom stats AIC
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel makeCluster
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#'
#' @rdname PermTestSurv
setMethod(f = "PermTestSurv", signature = "OmicsSurv",
          definition = function(OmicsSurv,
                                pathwayPCs_ls,
                                numReps = 0L,
                                parallel = FALSE,
                                numCores = NULL,
                                ...){

            ###  Function Setup  ###
            permute_SurvFit <- function(pathwayPCs_mat,
                                        resp_Surv,
                                        numReps_int = numReps){
              # browser()

              ###  True Model  ###
              pathwayPCs_mat <- as.matrix(pathwayPCs_mat)
              
              # We have an instance where all loadings and PC-values can be
              #   identically 0 (See Issue #69), so we add a catch for this:
              if(sum(abs(pathwayPCs_mat)) < .Machine$double.eps){
                return(1)
              }
              
              true_mod <- coxph(resp_Surv ~ pathwayPCs_mat)

              ###  p-Values  ###
              # Switch between real and permutation p-value
              if(numReps_int == 0){
                # Real score-based p-value

                trueMod_summ <- summary(true_mod)
                out_num <- unname(trueMod_summ$sctest["pvalue"])

              } else {
                # Permutation p-value

                permuteAIC_fun <- function(){

                  perm_resp <- SampleSurv(
                    response_vec = resp_Surv[, "time"],
                    event_vec = resp_Surv[, "status"],
                    parametric = FALSE
                  )

                  perm_Surv <- Surv(
                    time = perm_resp$response_vec,
                    event = perm_resp$event_vec
                  )

                  AIC(coxph(perm_Surv ~ pathwayPCs_mat))

                }

                trueAIC <- AIC(true_mod)
                permAIC <- replicate(n = numReps_int, expr = permuteAIC_fun())
                out_num <- mean(permAIC < trueAIC)

              }


              ###  Return  ###
              out_num

            }

            ###  Computation  ###
            # browser()
            response <- Surv(
              time = OmicsSurv@eventTime,
              event = OmicsSurv@eventObserved
            )

            if(parallel){
              # browser()

              ###  Parallel Computing Setup  ###
              message("Initializing Computing Cluster: ", appendLF = FALSE)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(
                deparse(quote(response)),
                deparse(quote(numReps))
              )
              clusterExport(
                cl = clust,
                varlist = clustVars_vec,
                envir = environment()
              )
              invisible(
                clusterEvalQ(cl = clust, library(pathwayPCA))
              )
              message("DONE")

              ###  Extract PCs  ###
              message("Extracting Pathway p-Values in Parallel: ",
                      appendLF = FALSE)
              pValues_vec <- parSapply(
                cl = clust,
                pathwayPCs_ls,
                permute_SurvFit,
                resp_Surv = response,
                numReps_int = numReps
              )
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway p-Values Serially")
              pValues_vec <- vapply(
                pathwayPCs_ls,
                permute_SurvFit,
                resp_Surv = response,
                numReps_int = numReps,
                FUN.VALUE = numeric(1)
              )
              message("DONE")

            }

            pValues_vec

          })
