context("SuperPCA_pVals")

data("colon_pathwayCollection")
data("colonSurv_df")

###  Survival  ###
test_that("Survival SuperPCA_pVals has correct class", {
  
  colon_OmicsSurv <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, 1:3],
    respType = "surv"
  )
  
  expect_s3_class(
    SuperPCA_pVals(
      object = colon_OmicsSurv,
      parallel = TRUE,
      numCores = 2,
      adjustpValues = TRUE,
      adjustment = "BH"
    ), class = "superpcOut"
  )
  
})


###  Regression  ###
test_that("Regression SuperPCA_pVals gives correct class", {

  colon_OmicsReg <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, 1:2],
    respType = "reg"
  )

  expect_s3_class(
    SuperPCA_pVals(
      object = colon_OmicsReg,
      parallel = TRUE,
      numCores = 2,
      adjustpValues = TRUE,
      adjustment = "BH"
    ), class = "superpcOut"
  )

})


###  Categorical  ###
test_that("Categorical SuperPCA_pVals gives correct class", {

  colon_OmicsCateg <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, c(1, 3)],
    respType = "categ"
  )

  expect_s3_class(
    SuperPCA_pVals(
      object = colon_OmicsCateg,
      parallel = TRUE,
      numCores = 2,
      adjustpValues = TRUE,
      adjustment = "BH"
    ), class = "superpcOut"
  )

})