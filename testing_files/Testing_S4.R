
######  Load Package and Data  ################################################

library(methods)
library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")



######  Set Validity  #########################################################

# MIGRATED TO R/validClass_Omics.R FILE


######  Load S4 Classes  ######################################################

# Pathway Extraction only
testOmicsPath <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection
)


# Survival (numeric and factor response)
testOmicsSurv <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[, 1:2],
  respType = "survival"
)
testOmicsSurv <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[-1, 1:2],
  respType = "survival"
)


# Regression (continuous response)
testOmicsReg <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df$OS_time,
  respType = "regression"
)
testOmicsReg <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df$OS_time[-1],
  respType = "regression"
)


# Classification (factor response)
testOmicsCateg <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df$OS_event,
  respType = "categorical"
)
testOmicsCateg <- create_Omics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df$OS_event[-1],
  respType = "categorical"
)


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
