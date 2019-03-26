context("LoadOntoPCs")



test_that("LoadOntoPCs gives correct errors", {
  
  gitHubPath_char <- "https://raw.githubusercontent.com/lizhongliu1996/pathwayPCAdata/master/"
  colon_aespcOut <- readRDS(
    url(paste0(gitHubPath_char, "colon_aespcOut.rds"))
  )
  
  expect_error(LoadOntoPCs(), 
               "")

  
})