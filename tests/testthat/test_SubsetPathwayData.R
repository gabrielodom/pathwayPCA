context("SubsetPathwayData")

data("colonSurv_df")
data("colon_pathwayCollection")

colon_OmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[, -c(2:3)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, 1:3],
  respType = "s"
)

test_that("SubsetPathwayData gives correct errors", {
  
    expect_error(
    SubsetPathwayData(
      colon_OmicsSurv,
      "haha"
    ), "Pathway not found."
  )
  
})

test_that("SubsetPathwayData returns df with correct classes", {
  
  expect_s3_class(
    SubsetPathwayData(
      colon_OmicsSurv,
      "KEGG_RETINOL_METABOLISM"
    ), c("tbl_df", "tbl", "data.frame")
  )
  
})