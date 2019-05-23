context("WhichPathways")


# The main functionality of this function is tested with the Contains function
data("colon_pathwayCollection")


test_that("WhichPathways() checks the class of x", {
  expect_error(
    WhichPathways(5, "MAP4K5"),
    "Object must be or extend class 'pathwayCollection'."
  )
})


test_that("WhichPathways() warns if a symbol is not found", {
  expect_warning(
    WhichPathways(colon_pathwayCollection, "AAAAA"),
    "No sets found to contain the requested symbols."
  )
})


colonSubset_PC <- WhichPathways(colon_pathwayCollection, "MAP4K5")

test_that("WhichPathways() returns a pathwayCollection", {
  expect_is(
    colonSubset_PC,
    c("pathwayCollection", "list")
  )
})

test_that("WhichPathways() does not add a 'description' field", {
  expect_true(
    is.null(colonSubset_PC$description)
  )
})



test_PC <- CreatePathwayCollection(
  sets_ls = list(set1 = LETTERS, set2 = letters),
  TERMS = c("UPPER_CASE", "lower_case"),
  description = c("louder", "quieter"),
  setType = "genes"
)
testSubset_PC <- WhichPathways(test_PC, "A")

test_that("WhichPathways() preserves the 'description' field", {
  expect_false(
    is.null(testSubset_PC$description)
  )
})

test_that("WhichPathways() preserves the set names", {
  expect_equal(
    names(test_PC)[1],
    names(testSubset_PC)[1]
  )
})