#' Test pathways with supervised PCA
#'
#' @description Given a supervised \code{OmicsPath} object (one of
#'    \code{OmicsSurv}, \code{OmicsReg}, or \code{OmicsCateg}), extract the
#'    first \eqn{k} principal components (PCs) from each pathway-subset of the
#'    -Omics assay design matrix, test their association with the response
#'    matrix, and return a data frame of the adjusted \eqn{p}-values for each
#'    pathway.
#'
#' @param object An object of superclass \code{OmicsPathway} with a response
#'    matrix or vector.
#' @param n.threshold The number of bins into which to split the feature scores
#'   in the fit object returned internally by the \code{\link{superpc.train}}
#'   function to the \code{\link{pathway_tScores}} and
#'   \code{\link{pathway_tControl}} functions. Defaults to 20. Smaller values
#'   may result in less accurate pathway \eqn{p}-values while larger values
#'   increase computation time.
#' @param numPCs The number of PCs to extract from each pathway. Defaults to 1.
#' @param parallel Should the computation be completed in parallel? Defaults to
#'   \code{FALSE}.
#' @param numCores If \code{parallel = TRUE}, how many cores should be used for
#'   computation? Defaults to \code{NULL}.
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
#' @return A data frame with columns:
#' \itemize{
#'   \item{\code{pathways} : }{The names of the pathways in the \code{Omics*}}
#'     object (given in \code{object@@trimPathwayCollection$pathways}.)
#'   \item{\code{setsize} : }{The number of genes in each of the original
#'     pathways (given in the \code{object@@trimPathwayCollection$setsize}
#'     object).}
#'   \item{\code{terms} : }{The pathway description, as given in the
#'     \code{object@@trimPathwayCollection$TERMS} object.}
#'   \item{\code{rawp} : }{The unadjusted \eqn{p}-values of each pathway.}
#'   \item{\code{...} : }{Additional columns as specified through the
#'     \code{adjustment} argument.}
#' }
#'
#' The data frame will be sorted in ascending order by the method specified
#'   first in the \code{adjustment} argument. If \code{adjustpValues = FALSE},
#'   then the data frame will be sorted by the raw \eqn{p}-values. If you have
#'   the suggested \code{tidyverse} package suite loaded, then this data frame
#'   will print as a \code{\link[tibble]{tibble}}. Otherwise, it will print as
#'   a data frame.
#'
#' @details This is a wrapper function for the \code{\link{pathway_tScores}},
#'   \code{\link{pathway_tControl}}, \code{\link{weibullMix_optimParams}},
#'   \code{\link{weibullMix_pValues}}, and \code{\link{adjust_and_sort}}
#'   functions.
#'
#'   Please see our Quickstart Guide for this package:
#'   \url{https://gabrielodom.github.io/pathwayPCA/articles/C1-Quickstart_Guide.html}
#'
#' @seealso \code{\link{create_OmicsPath}}; \code{\link{create_OmicsSurv}};
#'    \code{\link{create_OmicsReg}}; \code{\link{create_OmicsCateg}};
#'    \code{\link{pathway_tScores}}; \code{\link{pathway_tControl}};
#'    \code{\link{weibullMix_optimParams}}; \code{\link{weibullMix_pValues}};
#'    \code{\link{adjust_and_sort}}
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#' @include createClass_OmicsSurv.R
#' @include createClass_OmicsReg.R
#' @include createClass_OmicsCateg.R
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
#'   colonSurv_pVals_df <- superPCA_pVals(
#'     object = colon_OmicsSurv,
#'     parallel = TRUE,
#'     numCores = 16,
#'     adjustpValues = TRUE,
#'     adjustment = c("Hoch", "SidakSD")
#'   )
#' }
#'
#' @rdname superPCA_pVals
setGeneric("superPCA_pVals",
           function(object,
                    n.threshold = 20,
                    numPCs = 1,
                    parallel = FALSE,
                    numCores = NULL,
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
             standardGeneric("superPCA_pVals")
           }
)

#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'
#' @rdname superPCA_pVals
setMethod(f = "superPCA_pVals", signature = "OmicsPathway",
          definition = function(object,
                                n.threshold = 20,
                                numPCs = 1,
                                parallel = FALSE,
                                numCores = NULL,
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

            ###  Extract Information from S4 Object  ###
            geneArray_df <- t(object@assayData_df)
            obj_class <- class(object)

            switch (obj_class,
                    OmicsSurv = {

                      if(anyNA(object@eventTime) ||
                         anyNA(object@eventObserved)){
                        stop("\n Missing values are not permitted in the response for Supervised PCA.")
                      }

                      eventTime <- matrix(object@eventTime, ncol = 1)
                      eventObserved <- matrix(object@eventObserved, ncol = 1)
                      response_mat <- cbind(eventTime, eventObserved)
                      responseType <- "survival"

                    },
                    OmicsReg = {

                      if(anyNA(object@response)){
                        stop("\n Missing values are not permitted in the response for Supervised PCA.")
                      }

                      response_mat <- matrix(object@response, ncol = 1)
                      responseType <- "regression"

                    },
                    OmicsCateg = {

                      if(anyNA(object@response)){
                        stop("\n Missing values are not permitted in the response for Supervised PCA.")
                      }

                      # Matrices in R cannot be factors, so I'm going to try something
                      #   very hacky: a matrix in R is any atomic vector with a dim
                      #   attribute. This even works for ordered factors.
                      response_mat <- object@response
                      dim(response_mat) <- c(length(response_mat), 1)
                      responseType <- "categorical"

                    }
            )


            # responseType <- match.arg(responseType)
            if(adjustpValues){
              adjustment <- match.arg(adjustment, several.ok = TRUE)
            }

            pathwayGeneSets_ls <- object@trimPathwayCollection
            paths_ls <- pathwayGeneSets_ls$pathways
            minFeats <- attr(paths_ls, "minFeatures")
            if(parallel){

              ###  Parallel Computing Setup  ###
              message("Initializing Computing Cluster: ", appendLF = FALSE)
              clust <- makeCluster(numCores)
              clustVars_vec <- c(deparse(quote(paths_ls)),
                                 deparse(quote(geneArray_df)),
                                 deparse(quote(response_mat)))
              clusterExport(cl = clust,
                            varlist = clustVars_vec,
                            envir = environment())
              invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
              message("DONE")

              # browser()

              ###  Matrix of Student's t Scores and Controls  ###
              message("Calculating Pathway Test Statistics in Parallel: ",
                      appendLF = FALSE)
              tScores_ls <- parLapply(
                cl = clust,
                paths_ls,
                pathway_tScores,
                geneArray_df = geneArray_df,
                response_mat = response_mat,
                responseType = responseType,
                n.threshold = n.threshold,
                numPCs = numPCs,
                min.features = minFeats
              )
              message("DONE")

              message("Calculating Pathway Critical Values in Parallel: ",
                      appendLF = FALSE)
              tControl_ls <- parLapply(
                cl = clust,
                paths_ls,
                pathway_tControl,
                geneArray_df = geneArray_df,
                response_mat = response_mat,
                responseType = responseType,
                n.threshold = n.threshold,
                numPCs = numPCs,
                min.features = minFeats
              )
              message("DONE")
              stopCluster(clust)

            } else {

              ###  Matrix of Student's t Scores and Controls  ###
              message("Calculating Pathway Test Statistics Serially: ",
                      appendLF = FALSE)
              tScores_ls <- lapply(
                paths_ls,
                pathway_tScores,
                geneArray_df = geneArray_df,
                response_mat = response_mat,
                responseType = responseType,
                n.threshold = n.threshold,
                numPCs = numPCs,
                min.features = minFeats
              )
              message("DONE")

              message("Calculating Pathway Critical Values Serially: ",
                      appendLF = FALSE)
              tControl_ls <- lapply(
                paths_ls,
                pathway_tControl,
                geneArray_df = geneArray_df,
                response_mat = response_mat,
                responseType = responseType,
                n.threshold = n.threshold,
                numPCs = numPCs,
                min.features = minFeats
              )
              message("DONE")

            }

            # browser()


            ###  Maximum t-Score for each Gene Pathway  ###
            absMax <- function(vec){
              vec[which.max(abs(vec))]
            }
            # Lily told me that we only care about the t-scores from the first
            #   PC, so we will only extract the absolute maxima from the first
            #   row of the t-value matrices
            tScoreMax_num <- sapply(tScores_ls, function(ls) absMax(ls$tscor[1, ]) )
            tControlMax_num <- sapply(tControl_ls, function(x) absMax(x[1, ]) )


            ###  Calculate Raw Pathway p-Values  ###
            message("Calculating Pathway p-Values: ", appendLF = FALSE)
            genesPerPathway_vec <- unlist(pathwayGeneSets_ls$setsize)
            genesPerPathway_vec <- genesPerPathway_vec[names(paths_ls)]
            optParams_vec <- weibullMix_optimParams(
              max_tControl_vec = tControlMax_num,
              pathwaySize_vec = genesPerPathway_vec,
              ...
            )
            pvalues_vec <- weibullMix_pValues(
              tScore_vec = tScoreMax_num,
              pathwaySize_vec = genesPerPathway_vec,
              optimParams_vec = optParams_vec
            )
            message("DONE")

            if(adjustpValues){
              message("Adjusting p-Values and Sorting Pathway p-Value Data Frame: ",
                      appendLF = FALSE)
            } else {
              message("Sorting Pathway p-Value Data Frame: ", appendLF = FALSE)
            }

            out_df <- adjust_and_sort(pVals_vec = pvalues_vec,
                                      genesets_ls = pathwayGeneSets_ls,
                                      adjust = adjustpValues,
                                      proc_vec = adjustment,
                                      ...)
            message("DONE")


            ###  Wrangle PCA Output  ###
            sortedScores_ls <- tScores_ls[out_df$pathways]
            PCs_ls <- lapply(sortedScores_ls, `[[`, "PCs_mat")
            loadings_ls <- lapply(sortedScores_ls, `[[`, "loadings")

            ###  Return  ###
            out_ls <- list(
              pVals_df    = out_df,
              PCs_ls      = PCs_ls,
              loadings_ls = loadings_ls
            )

            class(out_ls) <- c("superpcOut", "pathwayPCA", "list")
            out_ls

          })
