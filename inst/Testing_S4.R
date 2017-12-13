
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
#  extract (but not subset) the score matrix [that's fixed now]. 9 min, 47 sec
#  for full extraction and subsetting.
pathway_AESPCs_1323_ls <- aespca2_ls
# devtools::use_data(pathway_AESPCs_1323_ls)



######  (TEST) Permutation Test Timing  #######################################
respIIIC <- ovarianFiltered_df$Tumor_Stage_Ovary_FIGO == "IIIC"
data("pathway_AESPCs_1323_ls")

logisticAIC <- function(response, predictor, permute = FALSE){

  if(permute){
    response <- sample(response)
  }

  glm(response ~ predictor, family = binomial(link = "logit"))$aic

}

# trueAIC <- logisticAIC(response = respIIIC,
#                        predictor = pathway_AESPCs_1323_ls[[1]])
#
# a <- Sys.time()
# distAIC2 <- replicate(10000, logisticAIC(response = respIIIC,
#                                      predictor = pathway_AESPCs_1323_ls[[1]],
#                                      permute = TRUE))
# Sys.time() - a # 18.55 sec for 10k; 3 min, 1 sec for 100k
#
# # The p-vaue (smaller AIC = better AIC)
# mean(trueAIC < distAIC)  # p = 0.34062
# mean(trueAIC < distAIC2) # p = 0.3379
# # I think those are close enough to warrant using 10k reps

a <- Sys.time()
pathwayAIC_p <- sapply(aespca2_ls[1:10], function(path){

  trueAIC <- logisticAIC(response = respIIIC, predictor = path)

  AIC_perm <- replicate(10000,
                        logisticAIC(response = respIIIC,
                                    predictor = path,
                                    permute = TRUE))

  # The p-vaue (smaller AIC = better AIC)
  mean(trueAIC < AIC_perm)

})
Sys.time() - a # 3 min, 6 sec for first 10 pathways;


library(parallel)
clus <- makeCluster(detectCores() - 2)
clusterExport(cl = clus, varlist = ls())
clusterEvalQ(cl = clus, library(pathwayPCA))
a <- Sys.time()
pathwayAIC_p <- parSapply(cl = clus,
                          pathway_AESPCs_1323_ls,
                          function(path){

  trueAIC <- logisticAIC(response = respIIIC, predictor = path)

  AIC_perm <- replicate(10000,
                        logisticAIC(response = respIIIC,
                                    predictor = path,
                                    permute = TRUE))

  # The p-vaue (smaller AIC = better AIC)
  mean(trueAIC < AIC_perm)

})
Sys.time() - a # 3 min, 5 sec for first 132 pathways (10%). 29 min, 21 sec for
#   the whole pathway set.

# Adjust the p-values
aespca_pValues_df <- data.frame(PathNames = names(pathwayAIC_p),
                                raw_pVal = pathwayAIC_p,
                                adj_pVal = p.adjust(pathwayAIC_p, "BH"))
# devtools::use_data(aespca_pValues_df)



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
# Leave this commented out, to see what pieces we still need to code
# source("inst/superpc.txt")

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
# for(i in 1:length(geneset$pathways)){ # 1 hour, 13 minutes for ~1,300 pathways
for(i in 1:5){

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
                       n.PCs = 1,
                       min.features = 2,
                       n.threshold = 20)

  tscore[i,] <- st.obj$tscor

}
Sys.time() - a   # 90 sec for first 100 pathways; 6 min 42 sec for 500 pathways
#   1 hr, 47 min for 7,949 pathways
# My hypothesis is that it's not the number of pathways, but the number of
#   pathways in the tail that have a lot of genes in them. For the set of 7,949
#   pathways, we cut off the pathways with more than 180 genes. We expected the
#   calculations to take ~8k / ~1.3k = 6x longer, but they look only about 1.5x
#   longer [1.776612 / (13 / 60 + 1) ~= 1.46]. The main reason is the huge
#   pathways with 45+ (95th percentile) genes in them. The more of these "whale"
#   pathways we have, the longer the computation will take.
# Further, I added the threshold.ignore component: the first half of the
#   threshold values yield the exact same t-scores for most pathways. To speed
#   up computation, we should be able to skip these. The first 500 took 5 min
#   45 sec, barely a minute faster. Probably not worth it. We'll change the
#   default to 0%.

plot(tscore[1, ], ylim = c(min(tscore), max(tscore)), type = "l", lwd = 3)
for(i in 2:500){
  lines(tscore[i, ], col = colours()[i], lwd = 3)
}
first10 <- sapply(1:500, function(row){
  var(tscore[row, 1:10])
})
last10 <- sapply(1:500, function(row){
  var(tscore[row, 11:20])
})
boxplot(first10, last10)
# We really don't see a thresholding effect until the last half of the threshold
#   values. We could add a "greedy" option to ignore the quantiles below 50%.


# library(plyr)
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
                       n.PCs = 1,      # number of rows of tscor
                       min.features = 2,
                       n.threshold = 20) # number of columns

  st.obj$tscor  # This is where we break things: we need to change from tall
  #   to wide data to take advantage of the clean plyr approach.
  # EDIT: transposing here makes the data tidy, but we want it wide.

}
wrapper1_fun(geneset$pathways[[1]])
tScores_mat <- sapply(geneset$pathways[1:100], wrapper1_fun)
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

# Transpose the matrix to return it to "tall" form
tScores_mat <- t(tScores_mat)


###  Control t-Scores  ###
# We know that the t-scores don't actually follow a t-distribution anymore:
plot(density(rt(150000, df = 177)))
lines(density(as.vector(tScores_mat)))
# lines(density(rt(150000, df = 3)))
tScores_sd <- sd(as.vector(tScores_mat))
tScores_mean <- mean(as.vector(tScores_mat))
lines(density(rnorm(150000, mean = tScores_mean, sd = tScores_sd)))
# THey do follow a Normal distribution pretty well though.
