context("SubsetPathwayData")

gitHubPath_char <-
  "https://raw.githubusercontent.com/lizhongliu1996/pathwayPCAdata/master/"
colon_aespcOut <- readRDS(
  url(paste0(gitHubPath_char, "colon_aespcOut.rds"))
)

test_that("SubsetPathwayData gives correct errors", {
  
  
  
  expect_error(
    SubsetPathwayData(), 
    ""
  )
  
  
  
  
  
  
  
})