
######  Initial Tests  ########################################################
# Migrated from inst/Testing_S4.R


library(methods)
library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")

# # ** Source the functions in the aes_pca.R file first **
# # Also, the aespca function should run on the unscaled original pathways. Thus,
# #   these results will be drastically different from the svd() results above.
# aespca(testRedPath[[1]], n = nrow(testRedPath[[1]]), d = 3, type = "predictor")
# aespca(testRedPath[[2]], n = nrow(testRedPath[[2]]), type = "predictor")
#
# a <- Sys.time()
# aespca_ls <- lapply(testRedPath, aespca, n = 58, d = 3, type = "predictor")
# Sys.time() - a # 15 min, 7 sec for 1323 screened pathways
# # This certainly warrants some parallel options
#
# ###  Parallel AES-PCA  ###
# library(parallel)
# clus <- makeCluster(detectCores() - 2)
# clusterExport(cl = clus, varlist = "testRedPath")
# clusterEvalQ(cl = clus, library(pathwayPCA))
# a <- Sys.time()
# aespca2_ls <- parLapply(cl = clus, testRedPath, function(x){
#   aespca(x, n = 58, d = 3, type = "predictor")$score
# })
# Sys.time() - a
# # 9 min, 41 sec for 1323 screened pathways over 18 cores, and we still have to
# #  extract (but not subset) the score matrix [that's fixed now]. 9 min, 47 sec
# #  for full extraction and subsetting.
# pathway_AESPCs_1323_ls <- aespca2_ls
# # devtools::use_data(pathway_AESPCs_1323_ls)



######  Test aespca() Alone  ##################################################
# I am testing this function now that the package is complete so that I can
#   better understand how the function works and what it is returning.
library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")

colon_OmicsPath <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection
)
# colonClean_OmicsPath <- expressedOmes(colon_OmicsPath)
testPathway_vec <- colon_OmicsPath@trimPathwayCollection$pathways[[1]]
testPathway_df <- colon_OmicsPath@assayData_df[, testPathway_vec]

aespca(X = testPathway_df, d = 2)



######  Full Walkthrough  #####################################################
# ovarian_OmicsPath <- create_OmicsPath(assayData_df = ovarianFiltered_df[, -(1:3)],
#                                       pathwayCollection_ls = aespca_Genesets_ls)
#
# ovarian_Omes <- expressedOmes(ovarian_OmicsPath)

# EXTRACTED TO extract_aesPCs_from_OmicsPath.R



# Test
# extract_aesPCs(object = colon_OmicsPath)
a <- Sys.time()
pcs_ls <- extract_aesPCs(colon_OmicsPath,
                         parallel = TRUE,
                         numCores = 15)
Sys.time() - a
# 3.688258 min to print to screen; 3.321165 min to save;
# 3.262908 min for load balancing (the upside to LB is that the brunt of the
# computation was done very quickly.)

# a <- Sys.time()
# pcs2_ls <- extract_aesPCs(ovarian_OmicsPath,
#                           parallel = TRUE,
#                           numCores = detectCores() - 2)
# Sys.time() - a # 3.229954 min; so, it works...
#
# ###  Create Ovarian OmicsCateg Object  ###
# tumour_fact <- as.factor(ovarianFiltered_df$Tumor_Stage_Ovary_FIGO)
# ovarian_OmicsCateg <- create_OmicsCateg(assayData_df = ovarianFiltered_df[, -(1:3)],
#                                         pathwayCollection_ls = aespca_Genesets_ls,
#                                         response_fact = tumour_fact)
#
# b <- Sys.time()
# ovarian_pVals <- permTest_OmicsCateg(OmicsCateg = ovarian_OmicsCateg,
#                                      pathwayPCs_ls = pcs2_ls,
#                                      numReps = 1000,
#                                      parallel = TRUE,
#                                      numCores = detectCores() - 2)
# Sys.time() - b # 7.418932 min for 1000 reps on laptop



######  Test AESPCA_pVals()  ##################################################


library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")

colon_OmicsSurv <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, 1:2],
  respType = "surv"
)

# AES-PCA
a0 <- Sys.time()
colon_aespcOut <- AESPCA_pVals(
  object = colon_OmicsSurv,
  numPCs = 1,
  numReps = 5000,
  parallel = TRUE,
  numCores = 15,
  adjustment = "BY"
)
Sys.time() - a0 # 34 seconds


# Regular PCA
a1 <- Sys.time()
colon_pcOut <- AESPCA_pVals(
  object = colon_OmicsSurv,
  numPCs = 1,
  numReps = 5000,
  parallel = TRUE,
  numCores = 15,
  asPCA = TRUE,
  adjustment = "BY"
)
Sys.time() - a1 # 34 seconds
