context("getPathpVals")

gitHubPath_char <-
  "https://raw.githubusercontent.com/lizhongliu1996/pathwayPCAdata/master/"
colon_aespcOut <- readRDS(
  url(paste0(gitHubPath_char, "colon_aespcOut.rds"))
)

test_that("getPathpVals returns df with correct classes", {
  
  expect_s3_class(
    getPathpVals(colon_aespcOut,
                 numPaths = 5), 
    class = c("tbl_df", "tbl", "data.frame")
  )
  
  expect_s3_class(
    getPathpVals(colon_aespcOut,
                 alpha = 0.01), 
    class = c("tbl_df", "tbl", "data.frame")
  )
  
})
