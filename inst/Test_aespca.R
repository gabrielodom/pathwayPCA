# Migrated from inst/Testing_S4.R


library(methods)
library(pathwayPCA)
load("data/ovarianFiltered_df.rda")
load("data/genesets_ls.rda")

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
