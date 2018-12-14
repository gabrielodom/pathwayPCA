
######  Load Package and Data  ################################################

library(methods)
library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")



######  Set Validity  #########################################################

# MIGRATED TO R/validClass_Omics.R FILE


######  Load S4 Classes  ######################################################

# Pathway Extraction only
testOmicsPath <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection
)


# Survival (numeric and factor response)
testOmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[, 1:3],
  respType = "survival"
)
testOmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[-1, 1:3],
  respType = "survival"
)


# Regression (continuous response)
testOmicsReg <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[, 1:2],
  respType = "regression"
)
testOmicsReg <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[-1, 1:2],
  respType = "regression"
)


# Classification (factor response)
testOmicsCateg <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[, c(1, 3)],
  respType = "categorical"
)
testOmicsCateg <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[-1, c(1, 3)],
  respType = "categorical"
)


######  Create Functions for these S4 Classes  ################################

# MIGRATED TO R/createOmics_TrimPathwayCollection.R FILE




######  PC Extraction with SVD  ###############################################

# Migrated to inst/Test_SVD_package_options.R


######  (TEST) PC Extraction with AES-PCA  ####################################

# Migrated to inst/Test_aespca.R



######  (TEST) Permutation Test Timing  #######################################

# Migrated to inst/Test_PermutationTest_Benchmark.R



######  (TEST) PC Extraction with Supervised PCA  #############################

# Migrated to inst/Test_supervisedPCA.R
