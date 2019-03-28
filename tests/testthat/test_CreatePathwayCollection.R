context("CreatePathwayCollection")

test_that("CreatePathwayCollection gives correct errors", {
  expect_error(CreatePathwayCollection(
    sets_ls = "pathway",
    TERMS = "five", 
    setType = "pathway", 
    setsize = 1
  ), 
  "The sets_ls object must be a list of set membership vectors.")
  
  expect_error(CreatePathwayCollection(
    sets_ls = list(one = 1, five = 5),
    TERMS = list(five = "five"), 
    setType = "pathway", 
    setsize = 1
  ), 
  "The TERMS object must be an atomic vector.")
  
  expect_error(CreatePathwayCollection(
    sets_ls = list(one = 1, five = 5),
    TERMS = "five", 
    setType = "gene", 
    setsize = 1
  ), 
  "The TERMS vector must have the same length as the sets list.")
})

test_that("CreatePathwayCollection gives correct warnings", {
  expect_warning(CreatePathwayCollection(
    sets_ls = list(one = 1,
                   two = 1:2,
                   three = 1:3,
                   four = 1:4,
                   five = 1:5),
    TERMS = c("one", "two", "three", "four", "five"),
    setType = "gene",
    setsize = 1
    ),
    "The names 'setsize' and 'n_tested' are reserved names of a pathwayCollection object.
            Values stored with these names may be overwritten during pathwayPCA function execution.
            Use with extreme caution.")
  
})
