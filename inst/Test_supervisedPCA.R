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
data("supervised_Genesets_ls")
geneset <- supervised_Genesets_ls


## run superpc test
# Leave this commented out, to see what pieces we still need to code
# source("inst/superpc.txt")

# This example is survival analysis, thus pinfo have both survival time and
#   censor status

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
# devtools::use_data(tScores_mat)


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

library(parallel)
clus <- makeCluster(detectCores() - 2)
clusterExport(cl = clus, varlist = ls())
clusterEvalQ(cl = clus, library(pathwayPCA))
a <- Sys.time()
tControl_mat <- parSapply(cl = clus, geneset$pathways, wrapperCtrl_fun)
Sys.time() - a # 9 min 3 sec

tControl_mat <- t(tControl_mat)
# Now we have the t-scores if the responses were random.
# devtools::use_data(tControl_mat)


######  Extreme Distribution and p-Values  ####################################

# Load the data files we need
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_Genesets_ls")
geneset <- supervised_Genesets_ls
data("tScores_mat")
data("tControl_mat")
rm(supervised_Tumors_df, supervised_Genesets_ls)


glen <- unlist(geneset$setsize)

# Find the largest and smallest t-scores
aa <- apply(tScores_mat, 1, max)
bb <- apply(tScores_mat, 1, min)

an1 <- sqrt(2 * log(glen))
top <- log(4 * pi) + log(log(glen))
bottom <- 2 * log(glen)
bn1 <- an1 * (1 - 0.5 * top / bottom)

newt <- aa
# If the negative number is larger in absolute value than the positive number,
#   then replace the positive value with the negative value
for(i in 1:length(aa)){
  if(abs(aa[i]) < abs(bb[i])){
    newt[i] <- bb[i]
  }
}

# There has seriously got to be an easier way to do this...


aa<-apply(tscore_control,1,max)
bb<-apply(tscore_control,1,min)


newc<-aa
for ( i in 1: length(aa) ) {
  if ( abs(aa[i])< abs(bb[i]) ) { newc[i]<-bb[i] }
}











mix.obj<-function(p,x,an, bn)
{
  z1<-(x-bn-p[2])*an*p[3]
  z2<-(x+bn+p[4])*an*p[5]
  e<-(p[1]*an*p[3])*exp(-z1-exp(-z1))+((1-p[1])*an*p[5])*exp(z2-exp(z2))
  if (any(e<=0)) Inf else -sum(log(e))
}

pp<-newc[newc>0]
nn<-newc[newc<0]

p0<-c(p=0.5,u1=1,s1=0.5 ,u2=1,s2=0.5)

lmix<-deriv(
  ~-log((p*an*s1)*exp(-((x-bn-u1)*an*s1)-exp(-(x-bn-u1)*an*s1))+((1-p)*an*s2)*exp(((x+bn+u2)*an*s2)-exp((x+bn+u2)*an*s2))),
  c("p","u1","s1","u2","s2"),
  function(x,an,bn,p,u1,s1,u2,s2) NULL)

mix.gr<-function(p,x,an,bn) {
  u1<-p[2]
  s1<-p[3]
  u2<-p[4]
  s2<-p[5]
  p<-p[1]
  colSums(attr(lmix(x,an,bn,p,u1,s1,u2,s2),"gradient"))}

aa<-optim(p0,mix.obj,mix.gr,x=newc,an=an1, bn=bn1, method="BFGS")

par<-aa$par



tt<-newt

newp<-rep(0,length(tt))

for ( i in 1:length(tt)) {
  tt1<-par[1]*(1-exp(-exp(-(tt[i]-bn1[i]-par[2])*an1[i]*par[3])))+(1-par[1])*(exp(-exp((tt[i]+bn1[i]+par[4])*an1[i]*par[5])))
  tt2<- 1-par[1]*(1-exp(-exp(-(tt[i]-bn1[i]-par[2])*an1[i]*par[3])))-(1-par[1])*(exp(-exp((tt[i]+bn1[i]+par[4])*an1[i]*par[5])))
  newp[i]<-min(tt1,tt2)
}

ntest<-data.frame(names(geneset$pathways),glen, newp)

bh<-mt.rawp2adjp(ntest$newp, "BH")
## by<-mt.rawp2adjp(all$p, "BY")
adjustedP<-bh$adjp[order(bh$index),]
all<-cbind(ntest, adjustedP)
all1<-cbind(all,unlist(geneset$TERMS))
all1<-all1[order(all1$BH,all1$rawp),]
all1<-all1[,-3]
names(all1)<-c("goterms"," setsize", "rawp", "FDR", "terms")


write.csv(all1, "results.csv")
