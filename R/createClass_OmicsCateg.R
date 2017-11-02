#' An S4 class for categorical responses within an "OmicsPathway" object
#'
#' @description This creates the "OmicsCateg" class which extends the
#'   "OmicsPathway" master class.
#'
#' @slot massSpec An $N x p$ data frame with named columns
#' @slot pathwaySet A list of known gene pathways
#' @slot response A factor vector of length $N$: the dependent variable of a
#'   generalized linear model
#'
#' @importFrom methods new
#'
#' @include createClass_OmicsPath.R validClass_Omics.R
#' @seealso \code{"\link[=OmicsPathway-class]{OmicsPathway}"}
#'
#' @export

create_OmicsCateg <- setClass("OmicsCateg",
                              slots = c(response = "factor"),
                              validity = valid_OmicsCateg,
                              contains = "OmicsPathway")
