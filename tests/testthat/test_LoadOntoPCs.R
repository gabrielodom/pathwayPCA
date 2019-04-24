context("LoadOntoPCs")

data("colonSurv_df")

gitHubPath_char <-
  "https://raw.githubusercontent.com/lizhongliu1996/pathwayPCAdata/master/"
colon_aespcOut <- readRDS(
  url(paste0(gitHubPath_char, "colon_aespcOut.rds"))
)


test_that("LoadOntoPCs returns df with correct classes", {
  
  expect_s3_class(
    LoadOntoPCs(
      design_df = colonSurv_df,
      loadings_ls = colon_aespcOut$loadings_ls
    ),
    c("tbl_df", "tbl", "data.frame")
  )
  
})
