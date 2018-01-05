# Migrated from inst/Testing_S4.R


library(methods)
library(pathwayPCA)
load("data/ovarianFiltered_df.rda")
load("data/genesets_ls.rda")


respIIIC <- ovarianFiltered_df$Tumor_Stage_Ovary_FIGO == "IIIC"
# data("pathway_AESPCs_1323_ls")
load("~/GitHub/pathwayPCA/inst/pathway_AESPCs_1323_ls.rda")

logisticAIC <- function(response, predictor, permute = FALSE){

  if(permute){
    response <- sample(response)
  }

  glm(response ~ predictor, family = binomial(link = "logit"))$aic

}

# trueAIC <- logisticAIC(response = respIIIC,
#                        predictor = pathway_AESPCs_1323_ls[[1]])
#
# a <- Sys.time()
# distAIC2 <- replicate(10000, logisticAIC(response = respIIIC,
#                                      predictor = pathway_AESPCs_1323_ls[[1]],
#                                      permute = TRUE))
# Sys.time() - a # 18.55 sec for 10k; 3 min, 1 sec for 100k
#
# # The p-vaue (smaller AIC = better AIC)
# mean(trueAIC < distAIC)  # p = 0.34062
# mean(trueAIC < distAIC2) # p = 0.3379
# # I think those are close enough to warrant using 10k reps

# a <- Sys.time()
# pathwayAIC_p <- sapply(aespca2_ls[1:10], function(path){
#
#   trueAIC <- logisticAIC(response = respIIIC, predictor = path)
#
#   AIC_perm <- replicate(10000,
#                         logisticAIC(response = respIIIC,
#                                     predictor = path,
#                                     permute = TRUE))
#
#   # The p-vaue (smaller AIC = better AIC)
#   mean(trueAIC < AIC_perm)
#
# })
# Sys.time() - a # 3 min, 6 sec for first 10 pathways;


library(parallel)
clus <- makeCluster(detectCores() - 2)
clusterExport(cl = clus, varlist = ls())
clusterEvalQ(cl = clus, library(pathwayPCA))
a <- Sys.time()
pathwayAIC_p <- parSapply(cl = clus,
                          pathway_AESPCs_1323_ls,
                          function(path){

                            trueAIC <- logisticAIC(response = respIIIC, predictor = path)

                            AIC_perm <- replicate(10000,
                                                  logisticAIC(response = respIIIC,
                                                              predictor = path,
                                                              permute = TRUE))

                            # The p-vaue (smaller AIC = better AIC)
                            mean(trueAIC < AIC_perm)

                          })
Sys.time() - a # 3 min, 5 sec for first 132 pathways (10%). 29 min, 28 sec for
#   the whole pathway set.

# Adjust the p-values
adjustedP <- adjustRaw_pVals(pathwayAIC_p, "BH")

# identical(names(pathway_AESPCs_1323_ls), names(pathwayAIC_p))
pathSize <- sapply(names(pathwayAIC_p), function(x){
  length(genesets_ls[x][[1]])
})

aespcaPathwaypVals_df <- data.frame(PathNames = names(pathwayAIC_p),
                                    setsize = pathSize)
rownames(aespcaPathwaypVals_df) <- NULL
all.equal(as.character(aespcaPathwaypVals_df$PathNames),
          names(pathway_AESPCs_1323_ls))

aespcaPathwaypVals_df <- cbind(aespcaPathwaypVals_df, adjustedP)
aespcaPathwaypVals_df <- aespcaPathwaypVals_df[order(aespcaPathwaypVals_df$BH,
                                                     aespcaPathwaypVals_df$rawp), ]
# devtools::use_data(aespcaPathwaypVals_df)
