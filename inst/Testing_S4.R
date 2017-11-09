
######  Load Package and Data  ################################################

library(methods)
library(pathwayPCA)
load("data/ovarianFiltered_df.rda")
load("data/genesets_ls.rda")



######  Set Validity  #########################################################
# # MIGRATED TO validClass_Omics.R FILE
# valid_OmicsSurv <- function(object){
#
#   nX <- nrow(object@massSpec)
#   nY <- length(object@eventTime)
#   nDelta <- length(object@eventObserved)
#
#   if(nY != nX){
#     return("Number of massSpec rows must match number of response times.")
#   } else if(nDelta != nX){
#     return("Number of massSpec rows must match number of response censoring indicators.")
#   } else {
#     return(TRUE)
#   }
#
# }
#
# valid_OmicsReg <- function(object){
#
#   nX <- nrow(object@massSpec)
#   nY <- length(object@response)
#
#   if(nY != nX){
#     return("Number of massSpec rows must match number of responses.")
#   } else {
#     return(TRUE)
#   }
#
# }
#
# valid_OmicsCateg <- function(object){
#
#   nX <- nrow(object@massSpec)
#   nY <- length(object@response)
#
#   if(nY != nX){
#     return("Number of massSpec rows must match number of responses.")
#   } else {
#     return(TRUE)
#   }
#
# }



######  Load S4 Classes  ######################################################

# Pathway Extraction only
testOmicsPath <- create_OmicsPath(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls)


# Survival (numeric and factor response)
Y_time <- rnorm(58, mean = 78, sd = 6)
Y_event <- sample(c(FALSE, TRUE), 58, replace = TRUE, prob = c(0.2, 1 - 0.2))
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls,
                                  eventTime = Y_time,
                                  eventObserved = Y_event)
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls,
                                  eventTime = Y_time[-1],
                                  eventObserved = Y_event)
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls,
                                  eventTime = Y_time,
                                  eventObserved = Y_event[-1])


# Regression (continuous response)
Y_reg <- rnorm(58)
testOmicsReg <- create_OmicsReg(massSpec = ovarianFiltered_df[, -(1:3)],
                                pathwaySet = genesets_ls,
                                response = Y_reg)
testOmicsReg <- create_OmicsReg(massSpec = ovarianFiltered_df[, -(1:3)],
                                pathwaySet = genesets_ls,
                                response = Y_reg[-1])


# Classification (factor response)
Y_class <- factor(sample(c("A", "B", "C"), 58, replace = TRUE))
testOmicsCateg <- create_OmicsCateg(massSpec = ovarianFiltered_df[, -(1:3)],
                                    pathwaySet = genesets_ls,
                                    response = Y_class)
testOmicsCateg <- create_OmicsCateg(massSpec = ovarianFiltered_df[, -(1:3)],
                                    pathwaySet = genesets_ls,
                                    response = Y_class[-1])


######  Create Functions for these S4 Classes  ################################

# ###  Geneset Extraction Function  ###
# # MIGRATED TO subsetExpressed-omes.R FILE
# setGeneric("expressedOmes",
#            function(object, ...){
#              standardGeneric("expressedOmes")
#              }
#            )
# setMethod(f = "expressedOmes", signature = "OmicsPathway",
#           definition = function(object, trim = 3){
#             genelist <- colnames(object@massSpec)
#
#             # Delete the genes from the pathways if they aren't recorded in our
#             #   data matrix
#             newPaths <- sapply(object@pathwaySet, function(x){
#               x[x %in% genelist]
#             })
#
#             # Remove any pathway that now has fewer than "trim" genes
#             newPaths_trim <- sapply(seq_along(newPaths), function(i){
#
#               shortPath <- length(newPaths[[i]]) < trim
#               if(shortPath){
#                 NULL
#               } else {
#                 newPaths[[i]]
#               }
#
#             })
#
#             names(newPaths_trim) <- names(object@pathwaySet)
#
#             nullPaths <- which(sapply(newPaths_trim, is.null))
#             paths_ls <- newPaths_trim[-nullPaths]
#
#             extractedMatrix_ls <- lapply(paths_ls, function(x){
#               object@massSpec[x]
#             })
#
#
#             attr(extractedMatrix_ls,
#                  "missingPaths") <- names(newPaths_trim)[nullPaths]
#             extractedMatrix_ls
#
#           })

# Test:
testRedPath <- expressedOmes(testOmicsPath)
testRedPath2 <- expressedOmes(testOmicsSurv)
testRedPath3 <- expressedOmes(testOmicsReg)
testRedPath4 <- expressedOmes(testOmicsCateg)
identical(testRedPath, testRedPath2)
identical(testRedPath, testRedPath3)
identical(testRedPath, testRedPath4)
# It works


######  PC Extraction with SVD  ###############################################
testRedPath <- expressedOmes(testOmicsPath)

# Notes: the symmetric = TRUE argument saves about 0.7 sec
gram <- function(X_mat){

  X_mat <- as.matrix(X_mat)
  t(X_mat) %*% X_mat

}

a <- Sys.time()
grams_ls <- lapply(testRedPath, gram) # 0.98 seconds
eigen_ls <- lapply(grams_ls, eigen, symmetric = TRUE)  # 2.67 seconds
Sys.time() - a # 3.43 sec

a <- Sys.time()
svd_ls <- lapply(testRedPath, svd)
Sys.time() - a # 1.53 sec; 37.3 Mb

library(rsvd)
a <- Sys.time()
svd_ls <- lapply(testRedPath, rsvd, k = 3)
Sys.time() - a # 1.95 sec, 4.2 Mb: subsets to 3 dimensions

# library(bootSVD)
# a <- Sys.time()
# svd_ls <- lapply(testRedPath, function(x){
#   fastSVD(as.matrix(x), nv = 3)
# })
# Sys.time() - a # 2.55 sec. Bad.

###  Let's go straight to the list of singular vectors  ###
redPathsScaled <- lapply(testRedPath, scale, center = TRUE, scale = TRUE)
# PCA
a <- Sys.time()
grams_ls <- lapply(redPathsScaled, gram)
eigen_ls <- lapply(covs_ls, function(x){

  svd_x <- eigen(x, symmetric = TRUE)
  svd_x$vectors[, 1:3]

})
Sys.time() - a # 3.45 sec


# SVD
a <- Sys.time()
svd1_ls <- lapply(redPathsScaled, function(x){

  svd_x <- svd(x)
  svd_x$v[, 1:3]

})
Sys.time() - a # 1.16 sec; 1.01 seconds average after restart (0.86 sec)


# RSVD
a <- Sys.time()
svd2_ls <- lapply(redPathsScaled, function(x){

  svd_x <- rsvd(x, k = 3)
  svd_x$v

})
Sys.time() - a # 1.45 sec

# I think that pretty well settles it: we're going with the base svd() function


######  (TEST) PC Extraction with AES-PCA  ####################################
# ** Source the functions in the aes_pca.R file first **
# Also, the aespca function should run on the unscaled original pathways. Thus,
#   these results will be drastically different from the svd() results above.
aespca(testRedPath[[1]], n = nrow(testRedPath[[1]]), d = 3, type = "predictor")
aespca(testRedPath[[2]], n = nrow(testRedPath[[2]]), type = "predictor")

a <- Sys.time()
aespca_ls <- lapply(testRedPath, aespca, n = 58, d = 3, type = "predictor")
Sys.time() - a # 15 min, 7 sec for 1323 screened pathways
# This certainly warrants some parallel options

###  Parallel AES-PCA  ###
library(parallel)
clus <- makeCluster(detectCores() - 2)
clusterExport(cl = clus, varlist = "testRedPath")
clusterEvalQ(cl = clus, library(pathwayPCA))
a <- Sys.time()
aespca2_ls <- parLapply(cl = clus, testRedPath, function(x){
  aespca(x, n = 58, d = 3, type = "predictor")$score
})
Sys.time() - a
# 9 min, 41 sec for 1323 screened pathways over 18 cores, and we still have to
#  extract (but not subset) the score matrix.


######  (TEST) PC Extraction with Supervised PCA  #############################
# setwd("C:/Users/gjo15/Dropbox (BBSR)/Ban and Odom - Bioconductor Package")
# # Pathway Load
# load("Original Chen Code/gene_expression/geneset.RData")  # [5,299] 4,266 pathways
# # load("Code/Supervised PCA/geneset_20171108.RData")      # [2,179] 7,949 pathways
# # Data Load
# # NOTE: data is p x n!
# load("Original Chen Code/gene_expression/array.RData")
# patInfo_df <- read.csv("Original Chen Code/gene_expression/pinfo.csv")

###  Load Data  ###
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_patInfo_df")
data("supervised_Genesets_ls")
geneset <- supervised_Genesets_ls


## run superpc test
source("inst/superpc.txt")

## This example is survival analysis, thus pinfo have both survival time annd censor status,

survY_df <- supervised_patInfo_df[, c("SurvivalTime", "disease_event")]
rm(supervised_Tumors_df, supervised_Genesets_ls, supervised_patInfo_df)

tscore <- array(0, dim = c(length(geneset$pathways), 20))
rownames(tscore) <- names(geneset$pathways)

###  The Basic Idea  ###
# Supervised PCA works like this:
#   1. Compute univariate model regression coefficients for each feature. That
#      is, given a candidate model y ~ f(x) + e, fit p independent models - one
#      for each gene in X. (For the pathway version, that's each gene in the
#      current pathway.)
#   2. Construct a reduced data matrix from genes / features whose univariate
#      model statistics exceed a certain threshold (estimated by CV).
#   3. Compute the first k PCs from this reduced data matrix.
#   4. Estimate a prediction model for y based on these first k PCs. (For the
#      pathway attribution exercise, check the significance of the pathway based
#      on the k PCs.)

a <- Sys.time()
for(i in 1:100){
  # for(i in 1:length(geneset$pathways)){ # 1 hour, 13 minutes for ~1,300 pathways
  # browser()

  genenames <- geneset$pathways[[i]]
  # pathway<-array[genename,]
  # y<-pinfo$SurvivalTime
  # censor<-pinfo$disease_event
  # pathwaysize=dim(pathway)[1]

  data <- list(x = array[genenames, ],
               y = survY_df$SurvivalTime,
               censoring.status = survY_df$disease_event,
               featurenames = genenames)
  ## if binary or continuous outcome
  ## data<-list(x=pathway,y=y, featurenames=genename)

  train <- superpc.train(data, type = "survival")
  ## if binary outcome: train<-superpc.train(data, type="binary")
  ## if continuous outcome: train<-superpc.train(data, type="continuous")

  st.obj <- superpc.st(fit = train,
                       data = data,
                       n.components = 1,
                       min.features = 2,
                       n.threshold = 20)

  tscore[i,] <- st.obj$tscor

}
Sys.time() - a   # 90 sec for first 100 pathways

library(plyr)
a <- Sys.time()
wrapper1_fun <- function(path){
  # browser()

  data <- list(x = array[path, ],
               y = survY_df$SurvivalTime,
               censoring.status = survY_df$disease_event,
               featurenames = path)

  train <- superpc.train(data, type = "survival")

  st.obj <- superpc.st(fit = train,
                       data = data,
                       n.components = 1,
                       min.features = 2,
                       n.threshold = 20)

  st.obj$tscor  # This is where we break things: we need to change from tall
  #   to wide data to take advantage of the clean plyr approach.
  # EDIT: transposing here makes the data tidy, but we want it wide.

}
tScores_ls <- sapply(geneset$pathways[1:100], wrapper1_fun)
Sys.time() - a   # 89.1 seconds, but we have a list instead of a matrix. But,
#   now we can run it in parallel! Using plyr::ldply returns a tall data frame
#   in 88.5 seconds, but we have to use plyr's parallel setup. Because the plyr
#   parallel options are difficult to set up, we changed to an sapply() call,
#   which returns a wide matrix with gene names as the column names. This takes
#   88.9 seconds to run, which is on par with plyr, but without the new package
#   call (it works smoothly with the clusterApply() syntax as well.)

###  Parallel Supervised PCA  ###
library(parallel)
clus <- makeCluster(detectCores() - 2)
clusterExport(cl = clus, varlist = ls())
clusterEvalQ(cl = clus, library(pathwayPCA))
a <- Sys.time()
tScores_mat <- parSapply(cl = clus, geneset$pathways, wrapper1_fun)
Sys.time() - a # 8 sec for parSapply, 7.71 sec for Load Balancing parSapplyLB
#   for first 100 pathways. 9 min, 15 seconds for all ~8k pathways with LB; 9
#   min, 12 seconds for all pathways without LB.
