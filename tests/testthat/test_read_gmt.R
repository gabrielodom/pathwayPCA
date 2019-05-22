context("read_gmt")


# Test read_gmt on connections
# Gabriel Odom
# 2019-05-22

# Lily wasn't able to run read_gmt on a connection, because I forgot to write
#   that code :/

# I've added a new argument to read_gmt: nChar. Because the file.info() function
#   does not work on connections, we need another way to get the number of
#   total characters to read from the file. I'm not sure how to check the number
#   of characters in a connection without opening the connection, so this will
#   probably need some more work in the future. For now, I default the number
#   of characters to twice the size of the MSigDB C5 pathway collection (about
#   6k pathways and nearly 5Mb in size, or nearly 5 million characters; there
#   isn't a 1-1 correspondence between character size and file size on disk, but
#   it's close).



test_that("read_gmt gives correct errors", {
  expect_error(read_gmt(
    12,
    setType = "region",
    description = FALSE
  ), 
  "'file' must be a character string or connection.")
})



GitHub_path <- "https://raw.githubusercontent.com/gabrielodom/pathwayPCA/master/"
con_path <- paste0(GitHub_path, "inst/extdata/wikipathways_human_symbol.gmt")

test_that("read_gmt opens a connection", {
  expect_is(read_gmt(
    url(con_path),
    description = TRUE
  ),
  c("pathwayCollection", "list"))
})



dataDir_path <- system.file(
  "extdata", package = "pathwayPCA", mustWork = TRUE
)

test_that("read_gmt reads a file", {
  expect_is(read_gmt(
    paste0(dataDir_path, "/wikipathways_human_symbol.gmt"),
    description = TRUE
  ),
  c("pathwayCollection", "list"))
})

