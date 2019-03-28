context("AESPCA_pVals")

test_that("AESPCA_pVals gives correct class results", {
  
  data("colonSurv_df")
  data("colon_pathwayCollection")
  colon_OmicsSurv <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, 1:3],
    respType = "surv"
  )
  
  #aespca
  expect_s3_class(AESPCA_pVals(
    object = colon_OmicsSurv,
    numPCs = 2,
    numReps = 10,
    parallel = TRUE,
    numCores = 4,
    adjustment = "BY"
  ), class = "aespcOut")
  
  # Regular PCA
  expect_s3_class(AESPCA_pVals(
    object = colon_OmicsSurv,
    numPCs = 1,
    numReps = 10,
    parallel = TRUE,
    numCores = 4,
    asPCA = TRUE,
    adjustment = "BY"
  ), class = "aespcOut")
  
  
})