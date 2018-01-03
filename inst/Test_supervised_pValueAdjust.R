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

# devtools::use_data(spcaPathwayPvals_df)


######  Graphics  #############################################################
