# Migrated from inst/Testing_S4.R

######  Load Data  ############################################################
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_patInfo_df")

data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls


## run superpc test
# Leave this commented out, to see what pieces we still need to code
# source("inst/superpc.txt")

# Note that the data is p x n
survY_df <- supervised_patInfo_df[, c("SurvivalTime", "disease_event")]
rm(supervised_Tumors_df, supervised_Genesets4240_ls, supervised_patInfo_df)


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



######  Parallel Supervised PCA  ##############################################

###  Make the wrapper a true function  ###
# Extracted to superPC_pathway_tScores.R

# Test
pathway_tScores(pathway_vec = geneset$pathways[[1]],
                geneArray_df = array,
                response_mat = survY_df$SurvivalTime,
                responseType = "regression")
# It works.


library(parallel)
clust <- makeCluster(detectCores() - 2)
clusterExport(cl = clust, varlist = ls())
clusterEvalQ(cl = clust, library(pathwayPCA))

a <- Sys.time()
tScores_mat <- parSapply(cl = clust,
                         geneset$pathways,
                         pathway_tScores,
                         geneArray_df = array,
                         response_mat = survY_df$SurvivalTime,
                         responseType = "regression")
Sys.time() - a # 1.467127 min for the 4,240 pathways with sizes in [5, 175]


# Transpose the matrix to return it to "tall" form
tScores_mat <- t(tScores_mat)
# tScoresReg4240_mat <- tScores_mat
# devtools::use_data(tScoresReg4240_mat)



######  Control t-Scores  #####################################################
# We know that the t-scores don't actually follow a t-distribution anymore:
plot(density(rt(150000, df = 177)))
lines(density(as.vector(tScores_mat)))
# This is appears to be a mixture distribution of two normals. This makes sense
#   because test statistics are often more positive or more negative


pathway_tControl(pathway_vec = geneset$pathways[[1]],
                 geneArray_df = array,
                 response_mat = survY_df$SurvivalTime,
                 responseType = "regression")

# TO DO FOR 20180103: WRITE THESE FUNCTIONS TO A FILE IN R/
#   THEN RUN THE PARALLEL CODE FOR THE CONTROL DATA AND SAVE THE RESULTS TO DATA/



######  Parallel Control t-Scores  ############################################
rm(list = ls())
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_patInfo_df")
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls
survY_df <- supervised_patInfo_df[, c("SurvivalTime", "disease_event")]
rm(supervised_Tumors_df, supervised_Genesets4240_ls, supervised_patInfo_df)


library(parallel)
clust <- makeCluster(detectCores() - 2)
clusterExport(cl = clust, varlist = ls())
clusterEvalQ(cl = clust, library(pathwayPCA))

a <- Sys.time()
tControl_mat <- parSapply(cl = clust,
                          geneset$pathways,
                          pathway_tControl,
                          geneArray_df = array,
                          response_mat = survY_df$SurvivalTime,
                          responseType = "regression")
Sys.time() - a # 1.213776 min

tControl_mat <- t(tControl_mat)
# tControlReg4240_mat <- tControl_mat
# Now we have the t-scores if the responses were random.
# devtools::use_data(tControlReg4240_mat)



######  Extreme Distribution  #################################################

# RESUME HERE 20180207
# Also, we'll need to document those data sets...

# EXTRACTED TO: superPC_optimeWeibullParams.R

# Load the data files we need
data("supervised_Genesets4240_ls")
pathwaylength_vec <- unlist(supervised_Genesets4240_ls$setsize)
rm(supervised_Genesets4240_ls)
data("tScoresReg4240_mat")
data("tControlReg4240_mat")


###  Find the largest and smallest t-scores  ###
absMax <- function(vec){
  vec[which.max(abs(vec))]
}
tScore_max <- apply(tScoresReg4240_mat, MARGIN = 1, FUN = absMax)
tControl_max <- apply(tControlReg4240_mat, MARGIN = 1, FUN = absMax)
rm(tControlReg4240_mat, tScoresReg4240_mat)


###  Find the Optimal Gumbel EV Mixture Parameters  ###
pOptim <- OptimGumbelMixParams(
  max_tControl_vec = tControl_max,
  pathwaySize_vec = pathwaylength_vec
)

pOptim



######  Calculate and Adjust p-Values  ########################################

# EXTRACTED TO: superPC_optimeWeibull_pValues.R

# Because pOptim was calculated using the tControl_max vector, that is how we
#   adjust our tScores_max vector by the control data set.
newp <- GumbelMixpValues(
  tScore_vec = tScore_max,
  pathwaySize_vec = pathwaylength_vec,
  optimParams_vec = pOptim
)


###  Adjust the p-Values  ###
# Load the data files we need
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls
rm(supervised_Genesets4240_ls)

# These are the p-values per pathway as returned by the newP_fun() function
ntest <- data.frame(goterms = names(geneset$pathways),
                    setsize = pathwaylength_vec,
                    rawp = newp)
rownames(ntest) <- NULL

bh <- ControlFDR(ntest$rawp, "BH")

ntest$FDR <- bh[, 2]
ntest$terms <- unlist(geneset$TERMS)
ntest <- ntest[order(ntest$FDR, ntest$rawp), ]
spcaRegPathwayPvals_df <- ntest

devtools::use_data(spcaRegPathwayPvals_df)
# write.csv(ntest, "results.csv")
