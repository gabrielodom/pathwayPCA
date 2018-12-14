# Test the computing speeds of different pca / svd routines.
# Gabriel Odom
# 20180601


# My candidates for "faster" pca / svd routines are
?corpcor::fast.svd # what we currently use
?irlba::irlba
?svd::trlan.svd; ?svd::trlan.eigen
?RSpectra::svds; ?RSpectra::eigs


# I believe that calculating the eigendecomposition of the pathway matrix may
#   be faster for pathways where p << n. I believe that calculating the singular
#   value decomposition of the pathway matrix will be faster for pathways where
#   p > n. I have no idea what will happen if p is slightly less than n. The
#   gains we get from eigendecomposition may be offset by calculating the
#   covariance matrix.

p <- 5
n <- 10


x_mat <- matrix(rnorm(n * p), ncol = p)

eigen(cor(x_mat))$vectors
svd(scale(x_mat))$v
all.equal(abs(eigen(cor(x_mat))$vectors),
          abs(svd(scale(x_mat))$v))
# Equal up to a sign.



######  SVD Tests  ############################################################

###  corpcor  ###
corpcor::fast.svd(scale(x_mat))$v
all.equal(corpcor::fast.svd(scale(x_mat))$v,
          svd(scale(x_mat))$v)


###  irlba  ###
irlba::irlba(scale(x_mat), nv = p - 1)$v
svd(scale(x_mat))$v
# These are fine, but I don't want to be warned for very small pathways. Also,
#   I can't test if these are exactly equal because I can't ask for all of the
#   singular vectors. This may not be a problem though. Let's see how many
#   pathways are very small:
library(pathwayPCA)
data("wikipwsHS_Entrez_pathwayCollection")
p_int <- lengths(wikipwsHS_Entrez_pathwayCollection$pathways)
summary(p_int)
plot(density(p_int))
quantile(p_int, 0.10)
quantile(p_int, 0.90)
# 10% of the pathways from wikipathways have six or fewer genes. 75% have 55 or
#   fewer. 90% have 101 genes or fewer. I think this tells us that we just need
#   to concentrate on the eigendecomposition. I'm going to keep testing the svd
#   though.


###  svd  ###
svd::ztrlan.svd(scale(x_mat))
# this is the complex decomposition. It would be nice if they had bloody
#   mentioned that minor detail in the help files.
svd::trlan.svd(scale(x_mat))
# This doesn't return the right singular vectors, and it gives lots of warnings


###  RSpectra  ###
RSpectra::svds(scale(x_mat), k = 5)$v
svd(scale(x_mat))$v
all.equal(RSpectra::svds(scale(x_mat), k = 5)$v,
          svd(scale(x_mat))$v)
# This is TRUE. The RSpectra code runs, but yells at us that we should be using
#   the base::svd() function instead. I think I'm ok with this warning, as it
#   will rarely trigger.



######  Eigendecomposition Tests  #############################################

###  corpcor  ###
# No pure eigendecomposition function. The same function "does both" (their
#   documentation mentions that fast.svd is best for "thin" or "fat" matrices.)
corpcor::fast.svd(scale(x_mat))$v
eigen(cor(x_mat))$vectors
all.equal(abs(corpcor::fast.svd(scale(x_mat))$v),
          abs(eigen(cor(x_mat))$vectors))
# Equal up to sign


###  irlba  ###
irlba::partial_eigen(cor(x_mat), n = p - 1)$`vectors`
eigen(cor(x_mat))$vectors
# These are fine, but I don't want to be warned for very small pathways. Also,
#   I can't test if these are exactly equal because I can't ask for all of the
#   singular vectors.


###  svd  ###
svd::trlan.eigen(cor(x_mat))
# This doesn't work


###  RSpectra  ###
RSpectra::eigs(cor(x_mat), k = 5)$vectors
eigen(cor(x_mat))$vectors
all.equal(abs(RSpectra::eigs(cor(x_mat), k = 5)$vectors),
          abs(eigen(cor(x_mat))$vectors))
# Equal up to sign



######  Testing Summary  ######################################################
# What we currently have is corpcor::fast.svd. If we switch, then the new
#   method must be *considerably* faster than this, otherwise it's not worth it
# Further, the svd:: package just doesn't work for small p, which is at least
#   10% of the pathways. Also, the irlba:: package gives pretty annoying
#   warnings for small p, and does not let us calculate all the PCs for the
#   smallest pathways. Not ok.



######  Benchmark Simulation  #################################################

###  Setup  ###
data("wikipwsHS_Entrez_pathwayCollection")
p_int <- lengths(wikipwsHS_Entrez_pathwayCollection$pathways)
# min(p_int)
p_int <- p_int + 4 # now the smallest pathways have length 5
n <- 100

set.seed(5)
testMatrices_ls <- lapply(p_int, function(P){
  matrix(rnorm(n * P), ncol = P)
})
# Now we have a "real" (enough) example list of data matrices. Let's see who
#   comes out on top


###  The Functions  ###
baseEigen <- function(inert){

  lapply(testMatrices_ls, function(mat){

    PCs <- eigen(cor(mat))$vectors
    mat %*% PCs[, 1:5]

  })

}
# Test
x1_ls <- baseEigen()

library(corpcor)
corpcorEigen <- function(inert){

  lapply(testMatrices_ls, function(mat){

    PCs <- fast.svd(scale(mat))$v
    mat %*% PCs[, 1:5]

  })

}
# Test
x2_ls <- corpcorEigen()

library(RSpectra)
spectraEigen <- function(inert){

  lapply(testMatrices_ls, function(mat){

    PCs <- eigs(cor(mat), k = 5)$vectors
    mat %*% PCs

  })

}
# Test
x3_ls <- spectraEigen() # 5 warnings
x3_ls <- suppressWarnings(spectraEigen())


###  Check the Functions  ###
sapply(seq_along(x1_ls), function(idx){
  all.equal(x1_ls[[idx]], x2_ls[[idx]])
}) # Base and corpcor equal up to a sign
sapply(seq_along(x1_ls), function(idx){
  all.equal(x1_ls[[idx]], x3_ls[[idx]])
}) # Base and RSpectra equal up to a sign


###  Replicate  ###
# The first time is for seed 1, second for 2, etc.
a1 <- Sys.time()
test1 <- replicate(100, x1_ls <- baseEigen())
Sys.time() - a1 # 1.569057, 1.567833, 1.512844, 1.500949, 1.611059 min
baseTimes_num <- c(1.569057, 1.567833, 1.512844, 1.500949, 1.611059)

a2 <- Sys.time()
test2 <- replicate(100, x2_ls <- corpcorEigen())
Sys.time() - a2 # 1.685411, 1.576553, 1.562411, 1.649198, 1.662961 min
corpcorTimes_num <- c(1.685411, 1.576553, 1.562411, 1.649198, 1.662961)

a3 <- Sys.time()
test3 <- replicate(100, x3_ls <- suppressWarnings(spectraEigen()))
Sys.time() - a3 # 58.55, 55.78702, 57.13908, 58.96591, 58.99293 sec
RSpectraTimes_num <- c(58.55, 55.78702, 57.13908, 58.96591, 58.99293) / 60

# It looks like the base eigendecomposition is better than the corpcor fast.svd
#   function. Brutal. However, the RSpectra functions are much faster.

library(ggplot2)
times_df <- data.frame(Package = c(rep("base", 5),
                                   rep("corpcor", 5),
                                   rep("RSpectra", 5)),
                       Time = c(baseTimes_num,
                                corpcorTimes_num,
                                RSpectraTimes_num))
ggplot(data = times_df) +
  aes(x = Package, y = Time) +
  scale_y_continuous(limits = c(0, 2)) +
  geom_boxplot()

summary(aov(Time ~ Package, data = times_df))
plot(aov(Time ~ Package, data = times_df))
