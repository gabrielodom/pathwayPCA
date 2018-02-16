# Extract a Wrapper Function from Test_supervisedPCA_regression.R


######  Load Data  ############################################################
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_patInfo_df")
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls

# Note that the data is p x n
survY_df <- supervised_patInfo_df[, c("SurvivalTime", "disease_event")]
rm(supervised_Tumors_df, supervised_Genesets4240_ls, supervised_patInfo_df)



######  Create Wrapper Function  ##############################################
calculate_pathway_pvalues <- function(pathwayGeneSets_ls,
                                      geneArray_df,
                                      response_mat,
                                      responseType = c("survival",
                                                       "regression",
                                                       "classification"),
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

  ###  Parallel Computing Setup  ###
  print("Initializing Cluster")
  require(parallel)
  clust <- makeCluster(numCores)
  paths_ls <- pathwayGeneSets_ls$pathways
  clustVars_vec <- c(deparse(quote(paths_ls)),
                     deparse(quote(geneArray_df)),
                     deparse(quote(response_mat)))
  clusterExport(cl = clust,
                varlist = clustVars_vec,
                envir = environment())
  invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
  print("DONE")


  ###  Matrix of Student's t Scores and Controls  ###
  print("Calculating Pathway Test Statistics")
  tScores_mat <- parSapply(cl = clust,
                           paths_ls,
                           pathway_tScores,
                           geneArray_df = geneArray_df,
                           response_mat = response_mat,
                           responseType = responseType)
  tScores_mat <- t(tScores_mat)
  print("DONE")

  print("Calculating Pathway Critical Values")
  tControl_mat <- parSapply(cl = clust,
                            paths_ls,
                            pathway_tControl,
                            geneArray_df = geneArray_df,
                            response_mat = response_mat,
                            responseType = responseType)
  tControl_mat <- t(tControl_mat)
  print("DONE")
  stopCluster(clust)


  ###  Maximum t-Score for each Gene Pathway  ###
  absMax <- function(vec){
    vec[which.max(abs(vec))]
  }
  tScoreMax_vec <- apply(tScores_mat, MARGIN = 1, FUN = absMax)
  tControlMax_vec <- apply(tControl_mat, MARGIN = 1, FUN = absMax)


  ###  Calculate Raw Pathway p-Values  ###
  print("Calculating Pathway p-Values")
  genesPerPathway_vec <- unlist(pathwayGeneSets_ls$setsize)
  optParams_vec <- weibullMix_optimParams(max_tControl_vec = tControlMax_vec,
                                   pathwaySize_vec = genesPerPathway_vec)
  pvalues_vec <- weibullMix_pValues(tScore_vec = tScoreMax_vec,
                                    pathwaySize_vec = genesPerPathway_vec,
                                    optimParams_vec = optParams_vec)

  out_df <- data.frame(goterms = names(pathwayGeneSets_ls$pathways),
                       terms = unlist(pathwayGeneSets_ls$TERMS),
                       setsize = genesPerPathway_vec,
                       rawp = pvalues_vec,
                       stringsAsFactors = FALSE)
  rownames(out_df) <- NULL
  print("DONE")


  ###  Adjust Pathway p-Values for Multiple Comparisons  ###
  print("Adjusting Pathway p-Values for Multiple Comparisons")
  adjusted_df <- data.frame(adjustRaw_pVals(out_df$rawp, adjustment))
  adjusted_df <- adjusted_df[, -1, drop = FALSE]
  out_df <- cbind(out_df, adjusted_df)

  out_df <- out_df[order(out_df[, adjustment[1]], out_df$rawp), ]
  print("DONE")


  ###  Return  ###
  out_df

}

###  Tests  ###
calculate_pathway_pvalues(pathwayGeneSets_ls = geneset,
                          geneArray_df = array,
                          response_mat = survY_df,
                          responseType = "survival",
                          parallel = TRUE,
                          numCores = detectCores() - 2,
                          adjustpValues = TRUE,
                          adjustment = "BH")
# It works

calculate_pathway_pvalues(pathwayGeneSets_ls = geneset,
                          geneArray_df = array,
                          response_mat = survY_df$SurvivalTime,
                          responseType = "regression",
                          parallel = TRUE,
                          numCores = detectCores() - 2,
                          adjustpValues = TRUE,
                          adjustment = "BH")
# It works


######  Wrapper Function v 2  #################################################
calculate_pathway_pvalues <- function(pathwayGeneSets_ls,
                                      geneArray_df,
                                      response_mat,
                                      responseType = c("survival",
                                                       "regression",
                                                       "classification"),
                                      n.threshold = 20,
                                      numPCs = 1,
                                      min.features = 5,
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
                                      alpha = 0.05,
                                      ...){
  # browser()

  responseType <- match.arg(responseType)
  if(adjustpValues){
    adjustment <- match.arg(adjustment, several.ok = TRUE)
  }

  if(parallel){

    ###  Parallel Computing Setup  ###
    message("Initializing Cluster")
    require(parallel)
    clust <- makeCluster(numCores)
    paths_ls <- pathwayGeneSets_ls$pathways
    clustVars_vec <- c(deparse(quote(paths_ls)),
                       deparse(quote(geneArray_df)),
                       deparse(quote(response_mat)))
    clusterExport(cl = clust,
                  varlist = clustVars_vec,
                  envir = environment())
    invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
    message("DONE")


    ###  Matrix of Student's t Scores and Controls  ###
    message("Calculating Pathway Test Statistics in Parallel")
    tScores_mat <- parSapply(cl = clust,
                             paths_ls,
                             pathway_tScores,
                             geneArray_df = geneArray_df,
                             response_mat = response_mat,
                             responseType = responseType,
                             n.threshold = n.threshold,
                             numPCs = numPCs,
                             min.features = min.features)
    tScores_mat <- t(tScores_mat)
    message("DONE")

    message("Calculating Pathway Critical Values in Parallel")
    tControl_mat <- parSapply(cl = clust,
                              paths_ls,
                              pathway_tControl,
                              geneArray_df = geneArray_df,
                              response_mat = response_mat,
                              responseType = responseType,
                              n.threshold = n.threshold,
                              numPCs = numPCs,
                              min.features = min.features)
    tControl_mat <- t(tControl_mat)
    message("DONE")
    stopCluster(clust)

  } else {

    ###  Matrix of Student's t Scores and Controls  ###
    message("Calculating Pathway Test Statistics Serially")
    tScores_mat <- sapply(paths_ls,
                          pathway_tScores,
                          geneArray_df = geneArray_df,
                          response_mat = response_mat,
                          responseType = responseType,
                          n.threshold = n.threshold,
                          numPCs = numPCs,
                          min.features = min.features)
    tScores_mat <- t(tScores_mat)
    message("DONE")

    message("Calculating Pathway Critical Values Serially")
    tControl_mat <- sapply(paths_ls,
                           pathway_tControl,
                           geneArray_df = geneArray_df,
                           response_mat = response_mat,
                           responseType = responseType,
                           n.threshold = n.threshold,
                           numPCs = numPCs,
                           min.features = min.features)
    tControl_mat <- t(tControl_mat)
    message("DONE")

  }




  ###  Maximum t-Score for each Gene Pathway  ###
  absMax <- function(vec){
    vec[which.max(abs(vec))]
  }
  tScoreMax_vec <- apply(tScores_mat, MARGIN = 1, FUN = absMax)
  tControlMax_vec <- apply(tControl_mat, MARGIN = 1, FUN = absMax)


  ###  Calculate Raw Pathway p-Values  ###
  message("Calculating Pathway p-Values")
  genesPerPathway_vec <- unlist(pathwayGeneSets_ls$setsize)
  optParams_vec <- weibullMix_optimParams(max_tControl_vec = tControlMax_vec,
                                          pathwaySize_vec = genesPerPathway_vec,
                                          ...)
  pvalues_vec <- weibullMix_pValues(tScore_vec = tScoreMax_vec,
                                    pathwaySize_vec = genesPerPathway_vec,
                                    optimParams_vec = optParams_vec)

  out_df <- data.frame(goterms = names(pathwayGeneSets_ls$pathways),
                       terms = unlist(pathwayGeneSets_ls$TERMS),
                       setsize = genesPerPathway_vec,
                       rawp = pvalues_vec,
                       stringsAsFactors = FALSE)
  rownames(out_df) <- NULL
  message("DONE")


  if(adjustpValues){

    ###  Adjust Pathway p-Values for Multiple Comparisons  ###
    message("Adjusting Pathway p-Values for Multiple Comparisons")
    adjusted_df <- data.frame(adjustRaw_pVals(out_df$rawp,
                                              adjustment,
                                              alpha = alpha))
    adjusted_df <- adjusted_df[, -1, drop = FALSE]
    out_df <- cbind(out_df, adjusted_df)

    out_df <- out_df[order(out_df[, adjustment[1]], out_df$rawp), ]
    message("DONE")

  } else {
    out_df <- out_df[order(out_df$rawp), ]
  }


  ###  Return  ###
  out_df

}

###  Tests  ###
a <- Sys.time()
survTest_df <- calculate_pathway_pvalues(pathwayGeneSets_ls = geneset,
                                         geneArray_df = array,
                                         response_mat = survY_df,
                                         responseType = "survival",
                                         parallel = TRUE,
                                         numCores = detectCores() - 2,
                                         adjustpValues = TRUE,
                                         adjustment = c("BH", "Hoch"))
Sys.time() - a # 10.19451 min
# It works

a <- Sys.time()
regTest_df <- calculate_pathway_pvalues(pathwayGeneSets_ls = geneset,
                                        geneArray_df = array,
                                        response_mat = survY_df$SurvivalTime,
                                        responseType = "regression",
                                        parallel = TRUE,
                                        numCores = detectCores() - 2,
                                        adjustpValues = TRUE,
                                        adjustment = c("BH", "Hoch"))
Sys.time() - a # 2.704036 min
# It works

# I need to test passing some arguments to the weibullMix_optimParams() function,
#   but this looks good overall. As a note, I don't want people passing any
#   arguments to the internal optim() function, because it is very sensitive.
a <- Sys.time()
regTest_df <- calculate_pathway_pvalues(pathwayGeneSets_ls = geneset,
                                        geneArray_df = array,
                                        response_mat = survY_df$SurvivalTime,
                                        responseType = "regression",
                                        parallel = TRUE,
                                        numCores = detectCores() - 2,
                                        adjustpValues = TRUE,
                                        adjustment = c("BH", "Hoch"),
                                        initialVals = c(p = 0.02,
                                                        mu1 = 1, s1 = 0.5,
                                                        mu2 = 1, s2 = 0.5))
Sys.time() - a # 2.704036 min
# Dots work


a <- Sys.time()
classifTest_df <- calculate_pathway_pvalues(pathwayGeneSets_ls = geneset,
                                            geneArray_df = array,
                                            response_mat = survY_df$disease_event,
                                            responseType = "classification",
                                            parallel = TRUE,
                                            numCores = detectCores() - 2,
                                            adjustpValues = TRUE,
                                            adjustment = c("BH", "SidakSS"))
Sys.time() - a # 3.700761 min



######  Wrapper Function v 3  #################################################
superPCA_pathway_pvals <- function(pathwayGeneSets_ls,
                                   geneArray_df,
                                   response_mat,
                                   responseType = c("survival",
                                                    "regression",
                                                    "classification"),
                                   n.threshold = 20,
                                   numPCs = 1,
                                   min.features = 5,
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
                                   alpha = 0.05,
                                   ...){
  # browser()

  responseType <- match.arg(responseType)
  if(adjustpValues){
    adjustment <- match.arg(adjustment, several.ok = TRUE)
  }

  if(parallel){

    ###  Parallel Computing Setup  ###
    message("Initializing Cluster")
    require(parallel)
    clust <- makeCluster(numCores)
    paths_ls <- pathwayGeneSets_ls$pathways
    clustVars_vec <- c(deparse(quote(paths_ls)),
                       deparse(quote(geneArray_df)),
                       deparse(quote(response_mat)))
    clusterExport(cl = clust,
                  varlist = clustVars_vec,
                  envir = environment())
    invisible(clusterEvalQ(cl = clust, library(pathwayPCA)))
    message("DONE")


    ###  Matrix of Student's t Scores and Controls  ###
    message("Calculating Pathway Test Statistics in Parallel")
    tScores_mat <- parSapply(cl = clust,
                             paths_ls,
                             pathway_tScores,
                             geneArray_df = geneArray_df,
                             response_mat = response_mat,
                             responseType = responseType,
                             n.threshold = n.threshold,
                             numPCs = numPCs,
                             min.features = min.features)
    tScores_mat <- t(tScores_mat)
    message("DONE")

    message("Calculating Pathway Critical Values in Parallel")
    tControl_mat <- parSapply(cl = clust,
                              paths_ls,
                              pathway_tControl,
                              geneArray_df = geneArray_df,
                              response_mat = response_mat,
                              responseType = responseType,
                              n.threshold = n.threshold,
                              numPCs = numPCs,
                              min.features = min.features)
    tControl_mat <- t(tControl_mat)
    message("DONE")
    stopCluster(clust)

  } else {

    ###  Matrix of Student's t Scores and Controls  ###
    message("Calculating Pathway Test Statistics Serially")
    tScores_mat <- sapply(paths_ls,
                          pathway_tScores,
                          geneArray_df = geneArray_df,
                          response_mat = response_mat,
                          responseType = responseType,
                          n.threshold = n.threshold,
                          numPCs = numPCs,
                          min.features = min.features)
    tScores_mat <- t(tScores_mat)
    message("DONE")

    message("Calculating Pathway Critical Values Serially")
    tControl_mat <- sapply(paths_ls,
                           pathway_tControl,
                           geneArray_df = geneArray_df,
                           response_mat = response_mat,
                           responseType = responseType,
                           n.threshold = n.threshold,
                           numPCs = numPCs,
                           min.features = min.features)
    tControl_mat <- t(tControl_mat)
    message("DONE")

  }




  ###  Maximum t-Score for each Gene Pathway  ###
  absMax <- function(vec){
    vec[which.max(abs(vec))]
  }
  tScoreMax_vec <- apply(tScores_mat, MARGIN = 1, FUN = absMax)
  tControlMax_vec <- apply(tControl_mat, MARGIN = 1, FUN = absMax)


  ###  Calculate Raw Pathway p-Values  ###
  message("Calculating Pathway p-Values")
  genesPerPathway_vec <- unlist(pathwayGeneSets_ls$setsize)
  optParams_vec <- weibullMix_optimParams(max_tControl_vec = tControlMax_vec,
                                          pathwaySize_vec = genesPerPathway_vec,
                                          ...)
  pvalues_vec <- weibullMix_pValues(tScore_vec = tScoreMax_vec,
                                    pathwaySize_vec = genesPerPathway_vec,
                                    optimParams_vec = optParams_vec)

  if(adjustpValues){
    message("Adjusting p-Values and Sorting Pathway p-Value Data Frame")
  } else {
    message("Sorting Pathway p-Value Data Frame")
  }

  out_df <- adjust_and_sort(pVals_vec = pvalues_vec,
                            genesets_ls = pathwayGeneSets_ls,
                            adjust = adjustpValues,
                            proc_vec = adjustment,
                            ...)
  message("DONE")

  ###  Return  ###
  out_df

}

###  Tests  ###
a <- Sys.time()
survTest_df <- superPCA_pathway_pvals(pathwayGeneSets_ls = geneset,
                                      geneArray_df = array,
                                      response_mat = survY_df,
                                      responseType = "survival",
                                      parallel = TRUE,
                                      numCores = detectCores() - 2,
                                      adjustpValues = TRUE,
                                      adjustment = c("BH", "Hoch"))
Sys.time() - a # 10.19451 min
# It works

a <- Sys.time()
regTest_df <- superPCA_pathway_pvals(pathwayGeneSets_ls = geneset,
                                     geneArray_df = array,
                                     response_mat = survY_df$SurvivalTime,
                                     responseType = "regression",
                                     parallel = TRUE,
                                     numCores = detectCores() - 2,
                                     adjustpValues = TRUE,
                                     adjustment = c("BH", "Hoch"))
Sys.time() - a # 2.704036 min
# It works

# I need to test passing some arguments to the weibullMix_optimParams() function,
#   but this looks good overall. As a note, I don't want people passing any
#   arguments to the internal optim() function, because it is very sensitive.
a <- Sys.time()
regTest_df <- superPCA_pathway_pvals(pathwayGeneSets_ls = geneset,
                                     geneArray_df = array,
                                     response_mat = survY_df$SurvivalTime,
                                     responseType = "regression",
                                     parallel = TRUE,
                                     numCores = detectCores() - 2,
                                     adjustpValues = TRUE,
                                     adjustment = c("BH", "Hoch"),
                                     initialVals = c(p = 0.02,
                                                     mu1 = 1, s1 = 0.5,
                                                     mu2 = 1, s2 = 0.5))
Sys.time() - a # 2.704036 min
# Dots work


a <- Sys.time()
classifTest_df <- superPCA_pathway_pvals(pathwayGeneSets_ls = geneset,
                                         geneArray_df = array,
                                         response_mat = survY_df$disease_event,
                                         responseType = "classification",
                                         parallel = TRUE,
                                         numCores = detectCores() - 2,
                                         adjustpValues = TRUE,
                                         adjustment = c("BH", "SidakSS"))
Sys.time() - a # 3.700761 min


######  Gene Present Check  ###################################################
mean(geneset$pathways[[1]] %in% rownames(array))
unique(sapply(geneset$pathways, function(i){
  mean(i %in% rownames(array))
}))
# All genes in the pathway list are present in the columns of the matrix. That
#   means we've cleaned the gene set list to only include genes we have actually
#   observed.



######  Wrapper Function v 4  #################################################
superPCA_pathway_pvals <- function(Omics_object,
                                   n.threshold = 20,
                                   numPCs = 1,
                                   min.features = 3,
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
                                   alpha = 0.05,
                                   ...){
  # browser()

  ###  Extract Information from S4 Object  ###
  geneArray_df <- t(Omics_object@massSpec)
  pathwayGeneSets_ls <- Omics_object@pathwaySet
  obj_class <- class(Omics_object)
  switch (obj_class,
          OmicsSurv = {

            eventTime <- matrix(Omics_object@eventTime, ncol = 1)
            eventObserved <- matrix(Omics_object@eventObserved, ncol = 1)
            response_mat <- cbind(eventTime, eventObserved)
            responseType <- "survival"

          },
          OmicsReg = {

            response_mat <- matrix(Omics_object@response, ncol = 1)
            responseType <- "regression"

          },
          OmicsCateg = {

            # Matrices in R cannot be factors, so I'm going to try something
            #   very hacky: a matrix in R is any atomic vector with a dim
            #   attribute. This even works for ordered factors.
            response_mat <- Omics_object@response
            dim(response_mat) <- c(length(response_mat), 1)
            responseType <- "classification"

          }
  )


  # responseType <- match.arg(responseType)
  if(adjustpValues){
    adjustment <- match.arg(adjustment, several.ok = TRUE)
  }

  if(parallel){

    ###  Parallel Computing Setup  ###
    message("Initializing Cluster")
    require(parallel)
    clust <- makeCluster(numCores)
    paths_ls <- pathwayGeneSets_ls$pathways
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
    message("Calculating Pathway Test Statistics in Parallel")
    tScores_mat <- parSapply(cl = clust,
                             paths_ls,
                             pathway_tScores,
                             geneArray_df = geneArray_df,
                             response_mat = response_mat,
                             responseType = responseType,
                             n.threshold = n.threshold,
                             numPCs = numPCs,
                             min.features = min.features)
    tScores_mat <- t(tScores_mat)
    message("DONE")

    message("Calculating Pathway Critical Values in Parallel")
    tControl_mat <- parSapply(cl = clust,
                              paths_ls,
                              pathway_tControl,
                              geneArray_df = geneArray_df,
                              response_mat = response_mat,
                              responseType = responseType,
                              n.threshold = n.threshold,
                              numPCs = numPCs,
                              min.features = min.features)
    tControl_mat <- t(tControl_mat)
    message("DONE")
    stopCluster(clust)

  } else {

    ###  Matrix of Student's t Scores and Controls  ###
    message("Calculating Pathway Test Statistics Serially")
    tScores_mat <- sapply(paths_ls,
                          pathway_tScores,
                          geneArray_df = geneArray_df,
                          response_mat = response_mat,
                          responseType = responseType,
                          n.threshold = n.threshold,
                          numPCs = numPCs,
                          min.features = min.features)
    tScores_mat <- t(tScores_mat)
    message("DONE")

    message("Calculating Pathway Critical Values Serially")
    tControl_mat <- sapply(paths_ls,
                           pathway_tControl,
                           geneArray_df = geneArray_df,
                           response_mat = response_mat,
                           responseType = responseType,
                           n.threshold = n.threshold,
                           numPCs = numPCs,
                           min.features = min.features)
    tControl_mat <- t(tControl_mat)
    message("DONE")

  }




  ###  Maximum t-Score for each Gene Pathway  ###
  absMax <- function(vec){
    vec[which.max(abs(vec))]
  }
  tScoreMax_vec <- apply(tScores_mat, MARGIN = 1, FUN = absMax)
  tControlMax_vec <- apply(tControl_mat, MARGIN = 1, FUN = absMax)


  ###  Calculate Raw Pathway p-Values  ###
  message("Calculating Pathway p-Values")
  genesPerPathway_vec <- unlist(pathwayGeneSets_ls$setsize)
  optParams_vec <- weibullMix_optimParams(max_tControl_vec = tControlMax_vec,
                                          pathwaySize_vec = genesPerPathway_vec,
                                          ...)
  pvalues_vec <- weibullMix_pValues(tScore_vec = tScoreMax_vec,
                                    pathwaySize_vec = genesPerPathway_vec,
                                    optimParams_vec = optParams_vec)

  if(adjustpValues){
    message("Adjusting p-Values and Sorting Pathway p-Value Data Frame")
  } else {
    message("Sorting Pathway p-Value Data Frame")
  }

  out_df <- adjust_and_sort(pVals_vec = pvalues_vec,
                            genesets_ls = pathwayGeneSets_ls,
                            adjust = adjustpValues,
                            proc_vec = adjustment,
                            ...)
  message("DONE")

  ###  Return  ###
  out_df

}


###  Guts Checking  ###
# pathway_tScores(pathway_vec = paths_ls[[2]],
#                 geneArray_df = geneArray_df,
#                 response_mat = response_mat,
#                 responseType = responseType)
#
# try_tScores <- lapply(paths_ls, function(p){
#
#   tryCatch(pathway_tScores(pathway_vec = p,
#                            geneArray_df = geneArray_df,
#                            response_mat = response_mat,
#                            responseType = responseType),
#            error = function(e) NULL)
#
# })
# which(sapply(try_tScores, is.null))
# # Pathways 1257 and 1310 are null.
# paths_ls[[1257]] # 4 genes
# paths_ls[[1310]] # 4 genes
# # Both of these pathways have enough genes so that the pass the `trim` argument
# #   in the `extract_aesPCs()` or `expressedOmes()` functions, but too few for
# #   the `min.features = 5` argument in the `pathway_tScores()` and
# #   `pathway_tControl` functions (which pass this argument on to the internal
# #   `superpc.st()` function).
# pathway_tScores(pathway_vec = paths_ls[[1257]],
#                 geneArray_df = geneArray_df,
#                 response_mat = response_mat,
#                 responseType = responseType)
# pathway_tScores(pathway_vec = paths_ls[[1310]],
#                 geneArray_df = geneArray_df,
#                 response_mat = response_mat,
#                 responseType = responseType)
#
# tryP_tScores <- parLapply(cl = clust,
#                           paths_ls,
#                           function(p){
#
#                             tryCatch(pathway_tScores(pathway_vec = p,
#                                                      geneArray_df = geneArray_df,
#                                                      response_mat = response_mat,
#                                                      responseType = responseType,
#                                                      n.threshold = n.threshold,
#                                                      numPCs = numPCs,
#                                                      min.features = min.features),
#                                      error = function(e) NULL)
#
#                           })
# which(sapply(tryP_tScores, is.null))



######  V 4 Tests  ############################################################

# FUNCTION EXTRACTED TO: superPC_wrapper.R

library(parallel)
library(pathwayPCA)

data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_patInfo_df")
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls

# Note that the data is p x n
survY_df <- supervised_patInfo_df[, c("SurvivalTime", "disease_event")]
rm(supervised_Tumors_df, supervised_Genesets4240_ls, supervised_patInfo_df)


# REQUIRES THE TUMOUR SURVIVAL DATA SET
###  Tests  ###

tumour_OmicsSurv <- create_OmicsSurv(massSpec_df = as.data.frame(t(array)),
                                     pathwaySet_ls = geneset,
                                     eventTime_vec = survY_df$SurvivalTime,
                                     eventObserved_vec = as.logical(survY_df$disease_event))
a <- Sys.time()
survTest_df <- superPCA_pVals(object = tumour_OmicsSurv,
                              parallel = TRUE,
                              numCores = detectCores() - 2,
                              adjustpValues = TRUE,
                              adjustment = c("BH", "SidakSS"))
Sys.time() - a # 1.715719 min
# It works


tumour_OmicsReg <- create_OmicsReg(massSpec_df = as.data.frame(t(array)),
                                   pathwaySet_ls = geneset,
                                   response_num = survY_df$SurvivalTime)
a <- Sys.time()
regTest_df <- superPCA_pVals(object = tumour_OmicsReg,
                             parallel = TRUE,
                             numCores = detectCores() - 2,
                             adjustpValues = TRUE,
                             adjustment = c("BH", "SidakSS"))
Sys.time() - a # 1.054811 min
# It works


tumour_OmicsCateg <- create_OmicsCateg(massSpec_df = as.data.frame(t(array)),
                                       pathwaySet_ls = geneset,
                                       response_fact = as.factor(survY_df$disease_event))
a <- Sys.time()
classifTest_df <- superPCA_pVals(object = tumour_OmicsCateg,
                                 parallel = TRUE,
                                 numCores = detectCores() - 2,
                                 adjustpValues = TRUE,
                                 adjustment = c("BH", "SidakSS"))
Sys.time() - a # 2.206281 min

###  Guts Checks  ###
# pathway_tScores(pathway_vec = paths_ls[[2]],
#                 geneArray_df = geneArray_df,
#                 response_mat = response_mat,
#                 responseType = responseType)


# REQUIRES THE VANDERBILT COLON CANCER DATA SET
# ###  Survival Test  ###
# colon2_OmicsSurv <- expressedOmes(colon_OmicsSurv)
# a <- Sys.time()
# survTest_df <- superPCA_pathway_pvals(Omics_object = colon2_OmicsSurv,
#                                       parallel = TRUE,
#                                       numCores = detectCores() - 2,
#                                       adjustpValues = TRUE,
#                                       adjustment = c("BH", "Hoch"))
# Sys.time() - a # 1.261178 min
#
#
# ###  Regression Test  ###
# colon2_OmicsReg <- expressedOmes(colon_OmicsReg)
# b <- Sys.time()
# regTest_df <- superPCA_pathway_pvals(Omics_object = colon2_OmicsReg,
#                                      parallel = TRUE,
#                                      numCores = detectCores() - 2,
#                                      adjustpValues = TRUE,
#                                      adjustment = c("BH", "Hoch"))
# Sys.time() - b # 1.074979 min
#
#
# ### Classification Test  ###
# colon2_OmicsCateg <- expressedOmes(colon_OmicsCateg)
# c <- Sys.time()
# classifTest_df <- superPCA_pathway_pvals(Omics_object = colon2_OmicsCateg,
#                                          parallel = TRUE,
#                                          numCores = detectCores() - 2,
#                                          adjustpValues = TRUE,
#                                          adjustment = c("BH", "SidakSS"))
# Sys.time() - c # 1.851152 min


# REQUIRES THE PANG OVARIAN CANCER DATA SET
