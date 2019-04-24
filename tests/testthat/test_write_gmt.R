context("write_gmt")

pathwayCollection1 <- list(
  pathways = list(
    c("C1orf27", "NR5A1", "BLOC1S4", "C4orf50"),
    c("TARS2", "DUSP5", "GPR88")
  ),
  TERMS = c("C-or-f_paths", "randomPath2"),
  description = c("these are", "totally made up")
)
class(pathwayCollection1) <- c("pathwayCollection", "list")

toy_pathwayCollection <- list(
  pathways = list(
    c("C1orf27", "NR5A1", "BLOC1S4", "C4orf50"),
    c("TARS2", "DUSP5", "GPR88"),
    c("TRX-CAT3-1", "LINC01333", "LINC01499", "LINC01046", "LINC01149")
  ),
  TERMS = c("C-or-f_paths", "randomPath2"),
  description = c("these are", "totally made up", "pathways")
)
class(toy_pathwayCollection) <- c("pathwayCollection", "list")


test_that("write_gmt gives correct errors", {
  
  expect_error(
    write_gmt(pathwayCollection1),
    "argument \"file\" is missing, with no default"
  )
  
  expect_error(
    write_gmt(),
    "argument \"pathwayCollection\" is missing, with no default"
  )
  
  expect_error(
    write_gmt(toy_pathwayCollection, file = "example_pathway.gmt"),
    "Number of sets should match number of TERMS."
  )
  
})



