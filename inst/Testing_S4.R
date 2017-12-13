
######  Load Package and Data  ################################################

library(methods)
library(pathwayPCA)
load("data/ovarianFiltered_df.rda")
load("data/genesets_ls.rda")



######  Set Validity  #########################################################

# MIGRATED TO R/validClass_Omics.R FILE


######  Load S4 Classes  ######################################################

# Pathway Extraction only
testOmicsPath <- create_OmicsPath(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls)


# Survival (numeric and factor response)
Y_time <- rnorm(58, mean = 78, sd = 6)
Y_event <- sample(c(FALSE, TRUE), 58, replace = TRUE, prob = c(0.2, 1 - 0.2))
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls,
                                  eventTime = Y_time,
                                  eventObserved = Y_event)
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls,
                                  eventTime = Y_time[-1],
                                  eventObserved = Y_event)
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = genesets_ls,
                                  eventTime = Y_time,
                                  eventObserved = Y_event[-1])


# Regression (continuous response)
Y_reg <- rnorm(58)
testOmicsReg <- create_OmicsReg(massSpec = ovarianFiltered_df[, -(1:3)],
                                pathwaySet = genesets_ls,
                                response = Y_reg)
testOmicsReg <- create_OmicsReg(massSpec = ovarianFiltered_df[, -(1:3)],
                                pathwaySet = genesets_ls,
                                response = Y_reg[-1])


# Classification (factor response)
Y_class <- factor(sample(c("A", "B", "C"), 58, replace = TRUE))
testOmicsCateg <- create_OmicsCateg(massSpec = ovarianFiltered_df[, -(1:3)],
                                    pathwaySet = genesets_ls,
                                    response = Y_class)
testOmicsCateg <- create_OmicsCateg(massSpec = ovarianFiltered_df[, -(1:3)],
                                    pathwaySet = genesets_ls,
                                    response = Y_class[-1])


######  Create Functions for these S4 Classes  ################################

# MIGRATED TO R/subsetExpressed-omes.R FILE




######  PC Extraction with SVD  ###############################################

# Migrated to inst/Test_SVD_package_options.R


######  (TEST) PC Extraction with AES-PCA  ####################################

# Migrated to inst/Test_aespca.R



######  (TEST) Permutation Test Timing  #######################################

# Migrated to inst/Test_PermutationTest_Benchmark.R



######  (TEST) PC Extraction with Supervised PCA  #############################

# Migrated to inst/Test_supervisedPCA.R
