# migrated from inst/Testing_S4.R


library(methods)
library(pathwayPCA)
load("data/ovarianFiltered_df.rda")
load("data/genesets_ls.rda")

# Test and Compare PC Extraction

testRedPath <- IntersectOmicsPwyCollct(testOmicsPath)

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
