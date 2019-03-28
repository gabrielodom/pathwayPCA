context("CreateOmics")

test_that("CreateOmics gives correct errors", {
  
  expect_error(CreateOmics(), 
  "argument \"assayData_df\" is missing, with no default")
  
  expect_error(CreateOmics(response = "12"), 
  "Response type required when a response is given.")
  
  expect_error(CreateOmics(respType = "s"), 
  "Response must be specified for type = survival")

})

test_that("CreateOmics gives correct messages", {
  
  data("colonSurv_df")
  data("colon_pathwayCollection")
  
  # Pathway
  expect_message(
    CreateOmics(
      assayData_df = colonSurv_df[, -(2:3)],
      pathwayCollection_ls = colon_pathwayCollection
    ),
    "Creating object of class OmicsPathway."
  )

  # Survival
  expect_message(
    CreateOmics(
      assayData_df = colonSurv_df[, -c(2:3)],
      pathwayCollection_ls = colon_pathwayCollection,
      response = colonSurv_df[, 1:3],
      respType = "s"
    ),
    "Creating object of class OmicsSurv."
  )

  # Regression
  expect_message(
    CreateOmics(
      assayData_df = colonSurv_df[, -c(2:3)],
      pathwayCollection_ls = colon_pathwayCollection,
      response = colonSurv_df[, 1:2],
      respType = "r"
    ),
    "Creating object of class OmicsReg."
  )

  # Categorical
  expect_message(
    CreateOmics(
      assayData_df = colonSurv_df[, -c(2:3)],
      pathwayCollection_ls = colon_pathwayCollection,
      response = colonSurv_df[, c(1,3)],
      respType = "c"
    ),
    "Creating object of class OmicsCateg."
  )

})


