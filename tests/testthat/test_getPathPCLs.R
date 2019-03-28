context("getPathPCLs")

gitHubPath_char <-
  "https://raw.githubusercontent.com/lizhongliu1996/pathwayPCAdata/master/"
colon_aespcOut <- readRDS(
  url(paste0(gitHubPath_char, "colon_aespcOut.rds"))
)

###  Object Class  ###
# class checked in AESPCA /SuperPCA tests

###  Errors  ###
test_that("getPathPCLs gives correct errors", {
  
  expect_error(
    getPathPCLs(
      colon_aespcOut,
      ""
    ),
    "Supplied pathway does not match any pathway in the supplied object."
  )
  
})


###  Output Structure  ###
out_ls <- getPathPCLs(
  colon_aespcOut,
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY"
)

test_that("Structure of object returned by getPathPCLs", {
  
  expect_named(
    out_ls, 
    c("PCs", "Loadings", "pathway", "term", "description")
  )
  
  expect_s3_class(
    out_ls$PCs,
    c("tbl_df", "tbl", "data.frame")
  )
  
  expect_s3_class(
    out_ls$Loadings,
    c("tbl_df", "tbl", "data.frame")
  )
  
})