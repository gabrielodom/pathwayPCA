context("SE2Tidy")

library(SummarizedExperiment)
data(airway, package = "airway")

airway_df <- SE2Tidy(airway)

test_that("SE2Tidy transforms data correctly", {

  expect_s3_class(airway_df, "data.frame")
  # expect_equal(dim(airway_df), c(8, 64112))
  expect_equal(dim(airway_df), c(8, 63687))

})
