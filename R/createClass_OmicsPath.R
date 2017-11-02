#' An S4 class for mass spectronomy data and gene pathway lists
#'
#' @slot massSpec An $N x p$ data frame with named columns
#' @slot pathwaySet A list of known gene pathways
#'
#' @importFrom methods new
#'
#' @export
setClass("OmicsPathway",
         slots = c(massSpec = "data.frame", pathwaySet = "list"))
setOldClass(c("tbl_df", "tbl", "data.frame"))

#' @export
create_OmicsPath <- function(massSpec_df, pathwaySet_ls){
  new("OmicsPathway",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls)
}