#' AES-PCA permutation test of categorical response for pathway PCs
#'
#' @description Given an \code{OmicsCateg} object and a list of pathway PCs from
#'   the \code{\link{ExtractAESPCs}} function, test if each pathway with
#'   features recorded in the bio-assay design matrix is significantly related
#'   to the categorical response.
#'
#' @param OmicsCateg A data object of class \code{OmicsCateg}, created by the
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
#'   from each pathway and an object of class \code{OmicsCateg}. This function
#'   will then calculate the AIC of a multivariate generalized linear model (via
#'   the \code{\link[stats]{glm}} function with a \code{\link[stats]{binomial}}
#'   error family) with the original observations as response and the pathway
#'   principal components as the predictor matrix.
#'
#'   Then, this function will create \code{numReps} permutations of the
#'   categorical response, fit models to each of these permuted responses
#'   (holding the path predictor matrix fixed), and calculate the AIC of each
#'   model. This function will return a named vector of permutation
#'   \eqn{p}-values, where the value for each pathway is the proportion of
#'   models for which the AIC of the permuted response model is less than the
#'   AIC of the original model.  Note that the AIC and log-likelihood are
#'   proportional because the number of parameters in each pathway is constant.
#'
#'   In future versions, this function will also be able to calculate permuted
#'   \eqn{p}-values for multinomial logistic regression and proportional odds
#'   logistic regression models, for n-ary and ordered categorical responses,
#'   respectively.
#'
#' @seealso \code{\link{CreateOmics}}; \code{\link{ExtractAESPCs}};
#'   \code{\link[stats]{glm}}; \code{\link[stats]{binomial}};
#'   \code{\link{SampleCateg}}
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsCateg.R
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
#'     response = colonSurv_df[, c(1,3)],
#'     respType = "categ"
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
#'   PermTestCateg(
#'     OmicsCateg = colon_Omics,
#'     pathwayPCs_ls = colonPCs_ls$PCs,
#'     parallel = TRUE,
#'     numCores = 2
#'   )
#' }
#'
#' @rdname PermTestCateg
setGeneric("PermTestCateg",
           function(OmicsCateg,
                    pathwayPCs_ls,
                    numReps = 0L,
                    parallel = FALSE,
                    numCores = NULL,
                    ...){
             standardGeneric("PermTestCateg")
           }
)

#' @importFrom stats AIC
#' @importFrom stats anova
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel makeCluster
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#'
#' @rdname PermTestCateg
setMethod(f = "PermTestCateg", signature = "OmicsCateg",
          definition = function(OmicsCateg,
                                pathwayPCs_ls,
                                numReps = 0L,
                                parallel = FALSE,
                                numCores = NULL,
                                ...){

            ###  Function Setup  ###
            # The distribution object will be passed in either indirectly or
            #   directly through the OmicsCateg object. I will probably put some
            #   sort of class / ordering check in the CreateOmicsCateg()
            #   function to add a slot or attribute for "binary", "n_ary" or
            #   "ordered" responses. Then, all I need is a switch here to define
            #   the "dist" object.
            # Now that I think about it, this should probably be three seperate
            #   functions. Nope, because then I'd have to have a switch around
            #   the apply statements. No thanks. I don't like how inneficient
            #   this is.
            permute_CategFit <- function(pathwayPCs_mat,
                                         response,
                                         nullMod,
                                         numReps_int = numReps){
              # browser()

              ###  True Model  ###
              pathwayPCs_mat <- as.matrix(pathwayPCs_mat)
              true_mod <- glm(response ~ pathwayPCs_mat, family = binomial)

              ###  p-Values  ###
              # Switch between real and permutation p-value
              if(numReps_int == 0){
                # Real score-based p-value

                mod_anova <- anova(nullMod, true_mod, test = "Chisq")
                out_num <- mod_anova[2, "Pr(>Chi)"]

              } else {
                # Permutation p-value

                permuteAIC_fun <- function(){

                  perm_resp <- SampleCateg(
                    response,
                    parametric = FALSE
                  )

                  AIC(glm(perm_resp ~ pathwayPCs_mat, family = binomial))

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
            response <- OmicsCateg@response
            null_mod <- glm(response ~ 1, family = binomial)

            if(parallel){
              # browser()

              ###  Parallel Computing Setup  ###
              message("Initializing Computing Cluster: ", appendLF = FALSE)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(
                deparse(quote(response)),
                deparse(quote(numReps)),
                deparse(quote(null_mod))
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
                permute_CategFit,
                response = response,
                numReps_int = numReps,
                nullMod = null_mod
              )
              stopCluster(clust)
              message("DONE")

            } else {

              message("Extracting Pathway p-Values Serially: ",
                      appendLF = FALSE)
              pValues_vec <- vapply(
                pathwayPCs_ls,
                permute_CategFit,
                response = response,
                numReps_int = numReps,
                nullMod = null_mod,
                FUN.VALUE = numeric(1)
              )
              message("DONE")

            }

            pValues_vec

          })
