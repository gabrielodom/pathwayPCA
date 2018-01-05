###  Extreme Distribution Adjusted p-Values  ###
# Migrated from Test_supervisedPCA.R



######  Data Load  ############################################################
data("supervised_Tumors_df")
array <- supervised_Tumors_df
data("supervised_Genesets4240_ls")
geneset <- supervised_Genesets4240_ls
data("tScores4240_mat")
data("tControl4240_mat")
rm(supervised_Tumors_df, supervised_Genesets4240_ls)



######  Largest and Smallest Absolute t-Scores  ###############################
matrixAbsMax <- function(mat, byRow = TRUE){

  absMax <- function(vec){
    vec[which.max(abs(vec))]
  }

  if(byRow){
    apply(mat, MARGIN = 1, FUN = absMax)
  } else {
    apply(mat, MARGIN = 2, FUN = absMax)
  }

}


tScore_max <- matrixAbsMax(tScores4240_mat)
tControl_max <- matrixAbsMax(tControl4240_mat)



######  Setup for Likelihood Optimization  ####################################
# Exported to superPC_optimWeibullParams.R

# Test
pOptim <- weibullMix_optimParams(max_tControl_vec = tControl_max,
                                 pathwaySize_vec = unlist(geneset$setsize))





######  Do Things with the t-Scores  ##########################################
# What things? No idea.
# Clean, vectorised, and functional code that does some mystery calculation.

# EXPORTED TO superPC_pathway_pValues

# Test
spcaPathwayPvals_df <- pathway_pValues(optimParams_vec = pOptim,
                                       max_tScores_vec = tScore_max,
                                       genelist_ls = geneset,
                                       FDRadjust = TRUE,
                                       multTestProc = "BH")
# # Reorder the rows by p-value
# spcaPathwayPvals_df[order(spcaPathwayPvals_df$BH, spcaPathwayPvals_df$rawp), ]

# devtools::use_data(spcaPathwayPvals_df)


######  Graphics  #############################################################
head(spcaPathwayPvals_df)
library(magrittr)
library(tidyverse)
spcaPathwayPvals_df %<>% mutate(pathScoreRaw = -log(rawp))
spcaPathwayPvals_df %<>% mutate(pathScoreFDR = -log(FDR))
spcaPathwayPvals_df$rank <- 1:nrow(spcaPathwayPvals_df)
head(spcaPathwayPvals_df)
tail(spcaPathwayPvals_df)

library(ggplot2)
###  p-Values  ###
ggplot(data = spcaPathwayPvals_df,
       aes(x = rank, y = pathScoreRaw)) +
  geom_bar(stat = "identity")
ggplot(data = spcaPathwayPvals_df,
       aes(x = rank, y = pathScoreFDR)) +
  geom_bar(stat = "identity")

###  p-Values by Pathway Size  ###
ggplot(data = spcaPathwayPvals_df,
       aes(x = setsize, y = pathScoreRaw)) +
  geom_point()
ggplot(data = spcaPathwayPvals_df,
       aes(x = setsize, y = pathScoreFDR)) +
  geom_point()
# No relationship - good.

###  Top p-Values  ###
ggplot(data = spcaPathwayPvals_df,
       aes(x = rank, y = pathScoreRaw, colour = FDR)) +
  geom_bar(stat = "identity")
ggplot(data = spcaPathwayPvals_df[1:500,],
       aes(x = rank, y = pathScoreRaw, colour = FDR)) +
  geom_bar(stat = "identity")
ggplot(data = spcaPathwayPvals_df[1:100,],
       aes(x = rank, y = pathScoreRaw, colour = FDR)) +
  geom_bar(stat = "identity")

# I don't think I'll be able to build great graphs until I have som results
#   thank are actually significant. I'll ask James about this.
