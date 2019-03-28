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
  
  # Pathway
  expect_message(CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
                             pathwayCollection_ls = gene_set_ls),
                             "Creating object of class OmicsPathway.")

  # Survival
  expect_error(CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
                           pathwayCollection_ls = gene_set_ls,
                           response = joinedExperiment_df[, 2:3],
                           respType = "s"),
                           "Creating object of class OmicsSurv.")

  # Regression
  expect_error(CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
                           pathwayCollection_ls = gene_set_ls,
                           response = joinedExperiment_df[, 2],
                           respType = "r"),
                           "Creating object of class OmicsReg.")

  # Categorical
  expect_error(CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
                           pathwayCollection_ls = gene_set_ls,
                           response = joinedExperiment_df[, 3],
                           respType = "c"),
                           "Creating object of class OmicsCateg.")

})


