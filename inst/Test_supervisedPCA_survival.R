# Migrated from inst/Testing_S4.R

# setwd("C:/Users/gjo15/Dropbox (BBSR)/Ban and Odom - Bioconductor Package")
# # Pathway Load
# load("Original Chen Code/gene_expression/geneset.RData")  # [5,299] 4,266 pathways
# # load("Code/Supervised PCA/geneset_20171108.RData")      # [2,179] 7,949 pathways
# # Data Load
# # NOTE: data is p x n!
# load("Original Chen Code/gene_expression/array.RData")
# patInfo_df <- read.csv("Original Chen Code/gene_expression/pinfo.csv")

######  Load Data  ############################################################
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_patInfo_df")
# data("supervised_Genesets_ls")
# geneset <- supervised_Genesets_ls

data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls


## run superpc test
# Leave this commented out, to see what pieces we still need to code
# source("inst/superpc.txt")

# This example is survival analysis, thus pinfo have both survival time and
#   censor status

survY_df <- supervised_patInfo_df[, c("SurvivalTime", "disease_event")]
rm(supervised_Tumors_df, supervised_Genesets4240_ls, supervised_patInfo_df)


###  Remove Small and Large Pathway Sets  ###
# which(geneset$setsize >= 300) # No large pathways
# max(geneset$setsize)
# which(geneset$setsize < 5) # ~3,700 small pathways
# smallpaths_ind <- which(geneset$setsize < 5)
# genesetReduced <- list(pathways = geneset$pathways[-smallpaths_ind],
#                        TERMS = geneset$TERMS[-smallpaths_ind],
#                        setsize = geneset$setsize[-smallpaths_ind])
# # Steven's original gene set list included setsize values in [175, 300], so his
# #   list has 4266 pathways. Ours has 4240.
# # supervised_Genesets4240_ls <- genesetReduced
# # devtools::use_data(supervised_Genesets4240_ls)


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

# tscore <- array(0, dim = c(length(geneset$pathways), 20))
# rownames(tscore) <- names(geneset$pathways)
#
# a <- Sys.time()
# # for(i in 1:length(geneset$pathways)){ # 1 hour, 13 minutes for ~1,300 pathways
# for(i in 1:500){
#
#   # browser()
#
#   genenames <- geneset$pathways[[i]]
#   # pathway<-array[genename,]
#   # y<-pinfo$SurvivalTime
#   # censor<-pinfo$disease_event
#   # pathwaysize=dim(pathway)[1]
#
#   data <- list(x = array[genenames, ],
#                y = survY_df$SurvivalTime,
#                censoring.status = survY_df$disease_event,
#                featurenames = genenames)
#   ## if binary or continuous outcome
#   ## data<-list(x=pathway,y=y, featurenames=genename)
#
#   train <- superpc.train(data, type = "survival")
#   ## if binary outcome: train<-superpc.train(data, type="binary")
#   ## if continuous outcome: train<-superpc.train(data, type="continuous")
#
#   st.obj <- superpc.st(fit = train,
#                        data = data,
#                        n.PCs = 1,
#                        min.features = 2,
#                        n.threshold = 20)
#
#   tscore[i,] <- st.obj$tscor
#
# }
# Sys.time() - a   # 90 sec for first 100 pathways; 6 min 32 sec for 500 pathways
# #   1 hr, 47 min for 7,949 pathways
# # My hypothesis is that it's not the number of pathways, but the number of
# #   pathways in the tail that have a lot of genes in them. For the set of 7,949
# #   pathways, we cut off the pathways with more than 180 genes. We expected the
# #   calculations to take ~8k / ~1.3k = 6x longer, but they look only about 1.5x
# #   longer [1.776612 / (13 / 60 + 1) ~= 1.46]. The main reason is the huge
# #   pathways with 45+ (95th percentile) genes in them. The more of these "whale"
# #   pathways we have, the longer the computation will take.
# # Further, I added the threshold.ignore component: the first half of the
# #   threshold values yield the exact same t-scores for most pathways. To speed
# #   up computation, we should be able to skip these. The first 500 took 5 min
# #   45 sec, barely a minute faster. Probably not worth it. We'll change the
# #   default to 0%.
#
# plot(tscore[1, ], ylim = c(min(tscore), max(tscore)), type = "l", lwd = 3)
# for(i in 2:500){
#   lines(tscore[i, ], col = colours()[i], lwd = 3)
# }
#
# # ASIDE: STOP USING sapply() OVER ROWS / COLUMNS. That's what apply() is for.
# # first10 <- sapply(1:500, function(row){
# #   var(tscore[row, 1:10])
# # })
# first <- apply(tscore[1:500, ], MARGIN = 1, function(x){
#   var(x[1:10])
# })
#
# # last10 <- sapply(1:500, function(row){
# #   var(tscore[row, 11:20])
# # })
# last <- apply(tscore[1:500, ], MARGIN = 1, function(x){
#   var(x[11:20])
# })
#
# # boxplot(first10, last10)
# boxplot(first, last)
# # We really don't see a thresholding effect until the last half of the threshold
# #   values. We could add a "greedy" option to ignore the quantiles below 50%.
#


######  Parallel Supervised PCA  ##############################################
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


###  Make the wrapper a true function  ###
# Extracted to superPC_pathway_tScores.R

# Test
pathway_tScores(pathway_vec = geneset$pathways[[1]],
                geneArray_df = array,
                response_mat = survY_df,
                responseType = "survival")
# It works.


library(parallel)
clust <- makeCluster(detectCores() - 2)
clusterExport(cl = clust, varlist = ls())
clusterEvalQ(cl = clust, library(pathwayPCA))
a <- Sys.time()
tScores_mat <- parSapply(cl = clust, geneset$pathways, pathway_tScores)
Sys.time() - a # 8 sec for parSapply, 7.71 sec for Load Balancing parSapplyLB
#   for first 100 pathways. 9 min, 15 seconds for all ~8k pathways with LB; 9
#   min, 12 seconds for all pathways without LB.

a <- Sys.time()
tScores_mat <- parSapply(cl = clust,
                         geneset$pathways,
                         pathway_tScores,
                         geneArray_df = array,
                         response_mat = survY_df)
Sys.time() - a # 4.926993 min for the 4,240 pathways with sizes in [5, 175]


# Transpose the matrix to return it to "tall" form
tScores_mat <- t(tScores_mat)
# devtools::use_data(tScores_mat)
tScores4240_mat <- tScores_mat
# devtools::use_data(tScores4240_mat)

# a <- Sys.time()
# tScores2_mat <- parSapply(cl = clust, genesetReduced$pathways, wrapper1_fun)
# Sys.time() - a # 4.929381 min
# tScores2_mat <- t(tScores2_mat)



######  Control t-Scores  #####################################################
# We know that the t-scores don't actually follow a t-distribution anymore:
plot(density(rt(150000, df = 177)))
lines(density(as.vector(tScores_mat)))
# lines(density(rt(150000, df = 3)))
tScores_sd <- sd(as.vector(tScores_mat))
tScores_mean <- mean(as.vector(tScores_mat))
lines(density(rnorm(150000, mean = tScores_mean, sd = tScores_sd)))
# They do follow a Normal distribution pretty well though.

###  Steven's Code for the Control Distribution  ###
# run control distribution, I generate Weibull distribution from emprical
#   data we have

library(survival)
wei <- survreg(Surv(SurvivalTime,
                    disease_event) ~ 1,
               data = survY_df,
               dist = "weibull")

pCensor <- mean(survY_df$disease_event == 1)
nObs <- nrow(survY_df)
parametric <- FALSE

tControl_mat <- matrix(NA, nrow = length(geneset$pathways), ncol = 20)
rownames(tControl_mat) <- names(geneset$pathways)


# for(i in 1:length(geneset$pathways)){                    # 1 hour, 7 minutes
for(i in 1:5){
  # browser()

  if(parametric){
    times_vec <- rweibull(nObs, shape = 1 / wei$scale, scale = exp(wei$coef))
  } else {
    times_vec <- sample(survY_df$SurvivalTime)
  }

  censor_ind <- runif(nObs) < pCensor

  for(m in 1:nObs){

    if(censor_ind[m]){
      times_vec[m] <- runif(1, min = 0, max = times_vec[m])
    }

  }

  genename <- geneset$pathways[[i]]
  data <- list(x = array[genename, ],
               y = times_vec,
               censoring.status = censor_ind,
               featurenames = genename)

  train <- superpc.train(data, type = "survival")

  st.obj <- superpc.st(fit = train,
                       data = data,
                       n.PCs = 1,
                       min.features = 2,
                       n.threshold = 20)

  tControl_mat[i, ] <- st.obj$tscor
}

######  Parallel Controlled t-Scores  #########################################
# Now that we have a working for() loop, let's extract the internal function
wrapperCtrl_fun <- function(path){
  # browser()

  if(parametric){
    times_vec <- rweibull(nObs, shape = 1 / wei$scale, scale = exp(wei$coef))
  } else {
    times_vec <- sample(survY_df$SurvivalTime)
  }

  censor_ind <- runif(nObs) < pCensor

  for(m in 1:nObs){

    if(censor_ind[m]){
      times_vec[m] <- runif(1, min = 0, max = times_vec[m])
    }

  }

  data <- list(x = array[path, ],
               y = times_vec,
               censoring.status = censor_ind,
               featurenames = path)

  train <- superpc.train(data, type = "survival")

  st.obj <- superpc.st(fit = train,
                       data = data,
                       n.PCs = 1,
                       min.features = 2,
                       n.threshold = 20)

  st.obj$tscor

}
wrapperCtrl_fun(geneset$pathways[[1]])
a <- Sys.time()
tControl_mat <- sapply(geneset$pathways[1:500], wrapperCtrl_fun)
Sys.time() - a # 1min 47 sec for first 100 pathways, 7 min 50 sec for first 500


###  Make the wrapper a true function  ###
## Extracted to superPC_pathway_tControls.R

# After looking at this, pathway_tControl needs to call three distinct functions
#   based on the type of the response. Each function will

###  Survival  ###
# Exported to superPC_permuteSamples.R

# Test
sample_Survivalresp(response_vec = survY_df$SurvivalTime,
                    censor_vec = survY_df$disease_event)

###  Regression  ###
# Exported to superPC_permuteSamples.R

# Test
sample_Regresp(response_vec = survY_df$SurvivalTime)


###  Classification  ###
# Exported to superPC_permuteSamples.R

# Test
sample_Classifresp(response_vec = survY_df$disease_event)
sample_Classifresp(response_vec = c(rep("A", 20), rep("B", 10), rep("C", 15)))
sample_Classifresp(response_vec = as.factor(survY_df$disease_event))


###  The Control Data  ###
# Exported to superPC_pathway_tControl.R

# Test
pathway_tControl(pathway_vec = geneset$pathways[[1]],
                 geneArray_df = array,
                 response_mat = survY_df,
                 responseType = "survival")

# TO DO FOR 20180103: WRITE THESE FUNCTIONS TO A FILE IN R/
#   THEN RUN THE PARALLEL CODE FOR THE CONTROL DATA AND SAVE THE RESULTS TO DATA/


library(parallel)
clust <- makeCluster(detectCores() - 2)
clusterExport(cl = clust, varlist = ls())
clusterEvalQ(cl = clust, library(pathwayPCA))
a <- Sys.time()
tControl_mat <- parSapply(cl = clust, geneset$pathways, wrapperCtrl_fun)
Sys.time() - a # 9 min 3 sec

a <- Sys.time()
tControl_mat <- parSapply(cl = clust,
                          geneset$pathways,
                          pathway_tControl,
                          geneArray_df = array,
                          response_mat = survY_df)
Sys.time() - a # 4.948044 min for reduced pathway set


tControl4240_mat <- t(tControl_mat)
# Now we have the t-scores if the responses were random.
# devtools::use_data(tControl4240_mat)

# a <- Sys.time()
# tControl2_mat <- parSapply(cl = clust, genesetReduced$pathways, wrapperCtrl_fun)
# Sys.time() - a # 4.824962 min
#
# tControl2_mat <- t(tControl2_mat)


######  Extreme Distribution and p-Values  ####################################

# RESUME HERE 20180103

# Load the data files we need
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls
data("tScores4240_mat")
data("tControl4240_mat")
rm(supervised_Tumors_df, supervised_Genesets4240_ls)


###  Find the largest and smallest t-scores  ###

# # Scores matrix
# aa <- apply(tScores_mat, 1, max)
# bb <- apply(tScores_mat, 1, min)
#
# newt <- aa
# # If the negative number is larger in absolute value than the positive number,
# #   then replace the positive value with the negative value
# for(i in 1:length(aa)){
#   if(abs(aa[i]) < abs(bb[i])){
#     newt[i] <- bb[i]
#   }
# }

# There has seriously got to be an easier way to do this...
# My version
absMax <- function(vec){
  vec[which.max(abs(vec))]
}
tScore_max <- apply(tScores4240_mat, MARGIN = 1, FUN = absMax)

# I've never really used the apply() function itself. Boy howdy!
# tScore_max <- sapply(1:nrow(tScores_mat), function(row){
#   absMax(tScores_mat[row, ])
# })
# names(tScore_max) <- rownames(tScores_mat)
# plot(tScore_max, ylim = c(-6, 5))


# # Control Matrix
# aa<-apply(tscore_control,1,max)
# bb<-apply(tscore_control,1,min)
#
#
# newc<-aa
# for ( i in 1: length(aa) ) {
#   if ( abs(aa[i])< abs(bb[i]) ) { newc[i]<-bb[i] }
# }

tControl_max <- apply(tControl4240_mat, MARGIN = 1, FUN = absMax)
# tControl_max <- sapply(1:nrow(tControl_mat), function(row){
#   absMax(tControl_mat[row, ])
# })
# names(tControl_max) <- rownames(tControl_mat)
# plot(tControl_max, ylim = c(-6, 5))

###  Use for Comparing against Supervised PCA Results  ###
# tScore2_max <- apply(tScores2_mat, MARGIN = 1, FUN = absMax)
# tControl2_max <- apply(tControl2_mat, MARGIN = 1, FUN = absMax)



###  Setup for Big Derivative Calculation  ###
# # I honestly have no idea what this is for...
# glen <- unlist(geneset$setsize)
# an1 <- sqrt(2 * log(glen))
# top <- log(4 * pi) + log(log(glen))
# bottom <- 2 * log(glen)
# bn1 <- an1 * (1 - 0.5 * top / bottom)

# Let's make this a function of the length vector
calc_anbn <- function(length_vec){

  logn <- log(length_vec)
  an1 <- sqrt(2 * logn)
  top <- log(4 * pi) + log(logn)
  bottom <- 2 * logn
  bn1 <- an1 * (1 - 0.5 * top / bottom)

  list(an = an1, bn = bn1)

}

# Test
pathwaylength_vec <- unlist(geneset$setsize)
abn_ls <- calc_anbn(pathwaylength_vec)
# all.equal(abn_ls$an, an1)
# all.equal(abn_ls$bn, bn1)
# # We're good, just remember to replace an1 and bn1 in later calculations

###  More Supervised PCA Comparison  ###
# pathway2length_vec <- unlist(genesetReduced$setsize)
# abn2_ls <- calc_anbn(pathway2length_vec)



###  Some Other Function  ###

# I have no idea where these numbers come from:
p0 <- c(p = 0.5, u1 = 1, s1 = 0.5, u2 = 1, s2 = 0.5)
# but they appear to be the arguments of the function we're trying to minimize.

# I don't know what this does either, but it's the "innard_formula" object
#   written as a function and coded slightly differently
# Ok, so this function is equation (4) Chen et al (2008). It's a density, and
#   therefore should NEVER be negative. The formulation here is different from
#   the paper in the following
#   1. They define zi = (t - mu_i) / sigma_i
#   2. They divide p and (1 - p) by sigma 1 and 2, respectively.
# Further, p is a mixing proportion between two Gumbel Extreme Value pdfs, not
#   a p-value. This makes more sense, but it does not explain how the optimum
#   mixing proportion is outise [0,1]. The values u and s are placeholders for
#   the mean and standard deviation, respectively, but they are not used as
#   they should be: the s values are used more as precision (being multiplied
#   instead of used as a divisor).
mix.obj <- function(p_vec, maxt_vec, an_vec, bn_vec){
  # browser()

  z1 <- (maxt_vec - bn_vec - p_vec["u1"]) * an_vec * p_vec["s1"]
  z2 <- (maxt_vec + bn_vec + p_vec["u2"]) * an_vec * p_vec["s2"]
  e  <- (p_vec["p"] * an_vec * p_vec["s1"]) *
          exp(-z1 - exp(-z1)) +
        ((1 - p_vec["p"]) * an_vec * p_vec["s2"]) *
          exp(z2 - exp(z2))
  # if(any(e <= 0)) Inf else -sum(log(e))
  ifelse(test = any(e <= 0), yes = 10 ^ 200, no = -sum(log(e)))
  # if(any(e <= 0)) e <- e + min(e[e != 0])
  # -sum(log(e))

}
# Test?
mix.obj(p_vec = p0,
        maxt_vec = tControl_max,
        an_vec = abn_ls$an,
        bn_vec = abn_ls$bn)
# I don't know what this means, but it runs. I think it's a liklihood value at
#   a certain point

###  More Supervised PCA Comparison  ###
# mix.obj(p_vec = p0,
#         maxt_vec = tControl2_max,
#         an_vec = abn2_ls$an,
#         bn_vec = abn2_ls$bn)
# # 17077.14
# # These values are the optimal values returned by the optim() function in the
# #   Supervised PCA directory. We get different values in ours.
p_supr <- c(p = 0.4731049,
            u1 = -0.4838454,
            s1 = 0.6289701,
            u2 = -0.4176837,
            s2 = 0.6277370)
mix.obj(p_vec = p_supr,
        maxt_vec = tControl_max,
        an_vec = abn_ls$an,
        bn_vec = abn_ls$bn)
# Something is off. This value is 7,627.002 in the other directory, but nearly
#   twice that here: 13,702.79.
# What is different? an and bn are identically the same. tControl_max is totally
#   different from newc. The issue is that while these numbers are different,
#   their histograms are almost the same, meaning that the values I get here are
#   well-within the realm of possibility. How does I still get a p-value outside
#   [0,1]?


# # There aren't any values identically 0, but I guess we don't know that couldn't
# #   happen.
# pp <- tControl_max[tControl_max > 0]
# nn <- tControl_max[tControl_max < 0]
# # just kidding, we never use these values at all


###  Take the gradient of some hella complex formula  ###
mix.gradient <- function(p_vec, maxt_vec, an_vec, bn_vec){
  # browser()

  # The deriv() syntax expects
  #   deriv(expr, namevec, function.arg, ...)
  # We are inputting some god-awful long formula for expr. The namevec argument
  #   should be a character vector of the variables we will derive expr with
  #   respect to; we input c("p","u1","s1","u2","s2"). The function.arg argument
  #   "must be specified and non-NULL", so I'm not sure what we're doing there,
  #   seeing as we supplied a NULL function.

  innard_formula <- as.formula(~ -log(
      (p * an * s1) *
        exp(-((x - bn - u1) * an * s1) - exp(-(x - bn - u1) * an * s1)) +
      ((1 - p) * an * s2) *
        exp(((x + bn + u2) * an * s2) - exp((x + bn + u2) * an * s2))
                                     )
                              )

  lmix_fun <- deriv(innard_formula,
                    c("p", "u1", "s1", "u2", "s2"),
                    function(x, an, bn, p, u1, s1, u2, s2) NULL)

  gradient <- lmix_fun(x = maxt_vec,
                       an = an_vec, bn = bn_vec,
                       p = p_vec["p"],
                       u1 = p_vec["u1"], s1 = p_vec["s1"],
                       u2 = p_vec["u2"], s2 = p_vec["s2"])
  # Extract the gradient matrix from the gradient vector, then find the column
  #   sums of that matrix
  colSums(attr(gradient,"gradient"))

}
# Test? I still don't get what the point is though
mix.gradient(p_vec = p0,
             maxt_vec = tControl_max,
             an_vec = abn_ls$an,
             bn_vec = abn_ls$bn)

### More Supervised PCA Comparison  ###
# mix.gradient(p_vec = p0,
#              maxt_vec = tControl2_max,
#              an_vec = abn2_ls$an,
#              bn_vec = abn2_ls$bn)
# #           p          u1          s1          u2          s2
# #    68.56648  9826.18810 30719.46106  8915.99482 27638.08031

# This is the parameter value that optimizes the mix.obj function
pOptim <- optim(par = p0,
                fn = mix.obj,
                gr = mix.gradient,
                maxt_vec = tControl_max,
                an_vec = abn_ls$an,
                bn_vec = abn_ls$bn,
                method = "BFGS")$par
# So this optimization routine *must* have a NAMED VECTOR OF PARAMETERS, not a
#   list. Gradient specification required for method "BFGS".
# ISSUE: this value doesn't match what the original code yields. The unedited
#   code yields par =  0.4731049 -0.4838454  0.6289701 -0.4176837  0.6277370.
#   We probably won't hit this exactly, but we absolutely must have p in [0, 1].
pOptim2 <- optim(par = p0,
                 fn = mix.obj,
                 gr = mix.gradient,
                 maxt_vec = tControl2_max,
                 an_vec = abn2_ls$an,
                 bn_vec = abn2_ls$bn,
                 method = "BFGS")$par
# After changing the pathway set to (nearly) match Steven's original results, I
#   get very close to the same output from the test values of mix.obj() and
#   mix.gradient(), but the optimal parameters are still waaaaay off for p.
# I've checked the inputs: an, bn, and p0 are the exact same, but tControl_max
#   and tControl2_max (for the reduced gene set list) are totally different
#   from newc and each other. That said, the values returned by mix.obj() and
#   mix.gradient() are really close.


# RESUME WORK HERE:
# Figure out why p-value calculation is different in the Supervised PCA directory

mix.obj(p_vec = pOptim,
        maxt_vec = tControl_max,
        an_vec = abn_ls$an,
        bn_vec = abn_ls$bn)
# [1] -64246.66
# WHAT? This value is huge - and negative? The other values we've seen so far
#   have been positive. What gives? Am I missing a negative sign somewhere?
# First, define the mix.obj function from the Test_supervisedPCA_pvalues.R file
mix.obj(p = pOptim,
        x = tControl_max,
        an = abn_ls$an,
        bn = abn_ls$bn)
# [1] -64246.66
# The results are identical, so I didn't miss a sign anywhere. That means it's
#   not my code that's different, but the data.
# Redefine mix.obj

# Try putting a constraint on the bounds of p
pOptim3 <- optim(par = p0,
                 fn = mix.obj,
                 gr = mix.gradient,
                 maxt_vec = tControl_max,
                 an_vec = abn_ls$an,
                 bn_vec = abn_ls$bn,
                 method = "L-BFGS-B",
                 lower = c(0, -Inf, 0, -Inf, 0),
                 upper = c(1, Inf, Inf, Inf, Inf))$par
# This errors because fn yields Inf values. It's a proper density, so I have no
#   idea why that ifelse() statement in the mix.obj function would ever proc at
#   all. Someone put it there for a reason though
pOptim3

# The new bounds work!!!! I'm getting results very similar to the optimised
#   values from the Supervised PCA directory:
# Ours:
# p         u1         s1         u2         s2
# 0.4821345 -0.4448527  0.7056764 -0.3878479  0.6742516
# Supervised PCA directory:
# p         u1         s1         u2         s2
# 0.4731049 -0.4838454  0.6289701 -0.4176837  0.6277370
# Praise the Lord!!
# Anyway, I had to change the ifelse() statement to return 10 ^ 200 instead of
#   Inf, because the "L-BFGS-B" optimisation routine (the one allowing you to
#   constrain the parameters) requires the function argument supplied to fn =
#   to be uniformally bounded.



# pOptim <- aa$par
# formerly known as "par". Change this below as well. These object names are
#   killing me


###  Do Things with the t-Scores  ###
# What things? No idea.

# # Original Code:
# # pOptim = par; tScore_max = tt; abn_ls = (an1, bn1)
# tt<-newt
#
# newp<-rep(0,length(tt))
#
# for ( i in 1:length(tt)) {
#   tt1<-par[1]*(1-exp(-exp(-(tt[i]-bn1[i]-par[2])*an1[i]*par[3])))+(1-par[1])*(exp(-exp((tt[i]+bn1[i]+par[4])*an1[i]*par[5])))
#   tt2<- 1-par[1]*(1-exp(-exp(-(tt[i]-bn1[i]-par[2])*an1[i]*par[3])))-(1-par[1])*(exp(-exp((tt[i]+bn1[i]+par[4])*an1[i]*par[5])))
#   newp[i]<-min(tt1,tt2)
# }


# # Cleaned-up code
# tt <- tScore_max
# newp <- rep(0, length(tt))
# names(newp) <- names(tScore_max)
#
# for(i in 1:length(tt)){
#
#   A <- 1 - exp(-exp(-(tt[i] - abn_ls$bn[i] - pOptim["u1"]) * abn_ls$an[i] * pOptim["s1"]))
#   B <- exp(-exp((tt[i] + abn_ls$bn[i] + pOptim["u2"]) * abn_ls$an[i] * pOptim["s2"]))
#
#   tt1 <- pOptim["p"] * A + (1 - pOptim["p"]) * B
#   tt2 <- 1 - tt1
#
#   newp[i] <- min(tt1, tt2)
#
# }
# It's the last piece which isn't vectorised:
# all.equal(newp, min(tt1, tt2))
# all.equal(newp, apply(cbind(tt1, tt2), MARGIN = 1, FUN = min))
# It works

# # Clean, vectorised, and functional code
newP_fun <- function(tScore_vec, optimParams_vec, abCounts_ls){
  # browser()

  an_s1 <- abCounts_ls$an * optimParams_vec["s1"]
  arg_A <- -(tScore_vec - abCounts_ls$bn - optimParams_vec["u1"]) * an_s1
  A <- 1 - exp(-exp(arg_A)) # in [0, 1]

  an_s2 <- abCounts_ls$an * optimParams_vec["s2"]
  arg_B <- (tScore_vec + abCounts_ls$bn + optimParams_vec["u2"]) * an_s2
  B <- exp(-exp(arg_B)) # also in [0, 1]

  # These values will be between 0 and 1 if the "p" value is in [0, 1]
  tt1 <- optimParams_vec["p"] * A + (1 - optimParams_vec["p"]) * B
  tt2 <- 1 - tt1

  apply(cbind(tt1, tt2), MARGIN = 1, FUN = min)

}
# I don't really know what this function does, but I know that the output matches
#   Steven's original code.

# Because pOptim was calculated using the tControl_max vector, that is how we
#   adjust our tScores_max vector by the control data set.
newp <- newP_fun(tScore_vec = tScore_max,
                 optimParams_vec = pOptim3,
                 abCounts_ls = abn_ls)


# These are the p-values per pathway as returned by the newP_fun() function
ntest <- data.frame(goterms = names(geneset$pathways),
                    setsize = pathwaylength_vec,
                    rawp = newp)
rownames(ntest) <- NULL

bh <- multtest::mt.rawp2adjp(ntest$rawp, "BH")
adjustedP <- bh$adjp[order(bh$index), ]

ntest$FDR <- adjustedP[, 2]
ntest$terms <- unlist(geneset$TERMS)
ntest <- ntest[order(ntest$FDR, ntest$rawp), ]
spcaPathwayPvals_df <- ntest

devtools::use_data(spcaPathwayPvals_df)
# write.csv(ntest, "results.csv")
