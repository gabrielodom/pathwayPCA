context("Contains")


test_that("Contains() 'any' operates correctly.", {
  
  expect_true(
    Contains(
      long = LETTERS,
      short = c("A", "!"),
      matches = "any"
    )
  )
  
  expect_true(
    Contains(
      long = LETTERS,
      short = c("A", "B"),
      matches = "any"
    )
  )
  
})


test_that("Contains() defaults to 'any' match.", {
  expect_true(
    Contains(
      long = LETTERS,
      short = c("A", "!")
    )
  )
})


test_that("Contains() 'all' operates correctly.", {
  
  expect_false(
    Contains(
      long = LETTERS,
      short = c("A", "!"),
      matches = "all"
    )
  )
  
  expect_true(
    Contains(
      long = LETTERS,
      short = c("A", "B"),
      matches = "all"
    )
  )
  
})


test_that("Contains() partial matching works.", {
  
  expect_error(
    Contains(
      long = 1:10,
      short = c("A", "!"),
      matches = "all",
      partial = TRUE
    ),
    "Partial matching permitted only for character vectors."
  )
  
  expect_error(
    Contains(
      long = LETTERS,
      short = c("A", "!"),
      matches = "all",
      partial = TRUE
    ),
    "Partial matching permitted only for single values of 'short'."
  )
  
  expect_warning(
    Contains(
      long = LETTERS,
      short = "A",
      matches = "all",
      partial = TRUE
    ),
    "Partial matching takes any matches. matches = 'all' will be ignored."
  )
  
  expect_true(
    Contains(
      long = LETTERS,
      short = "A",
      matches = "any",
      partial = TRUE
    )
  )
  
})