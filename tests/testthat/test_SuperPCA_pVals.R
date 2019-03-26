context("SuperPCA_pVals")

test_that("SuperPCA_pVals gives correct class results", {
  
  data("colon_pathwayCollection")
  data("colonSurv_df")
  
  ######  Survival  #############################################################
  colon_OmicsSurv <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, 1:3],
    respType = "surv"
  )
  
  expect_s4_class(SuperPCA_pVals(
                  object = colon_OmicsSurv,
                  parallel = TRUE,
                  numCores = 15,
                  adjustpValues = TRUE,
                  adjustment = c("BH", "SidakSS")
                 ), class = "superpcOut")
  
  ######  Regression  ###########################################################
  colon_OmicsReg <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, 1:2],
    respType = "reg"
  )
  
  expect_s4_class(SuperPCA_pVals(
                  object = colon_OmicsReg,
                  parallel = TRUE,
                  numCores = 15,
                  adjustpValues = TRUE,
                  adjustment = c("BH", "SidakSS")
                 ), class = "superpcOut")
  
  ######  Categorical  ##########################################################
  colon_OmicsCateg <- CreateOmics(
    assayData_df = colonSurv_df[, -(2:3)],
    pathwayCollection_ls = colon_pathwayCollection,
    response = colonSurv_df[, c(1, 3)],
    respType = "categ"
  )
  
  expect_s4_class(SuperPCA_pVals(
                  object = colon_OmicsReg,
                  parallel = TRUE,
                  numCores = 15,
                  adjustpValues = TRUE,
                  adjustment = c("BH", "SidakSS")
                ), class = "superpcOut")
  
})