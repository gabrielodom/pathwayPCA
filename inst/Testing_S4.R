
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

a <- Sys.time()
# covs_ls <- lapply(testRedPath, cov) # 0.98 seconds
eigen_ls <- lapply(covs_ls, eigen, symmetric = TRUE)  # 2.67 seconds
Sys.time() - a # 3.43 sec

a <- Sys.time()
svd_ls <- lapply(testRedPath, svd)
Sys.time() - a # 1.53 sec
