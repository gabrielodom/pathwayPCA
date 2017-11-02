#' An S4 class for continuous responses within an "OmicsPathway" object
#'
#' @description This creates the "OmicsReg" class which extends the
#'   "OmicsPathway" master class.
#'
#' @slot massSpec An $N x p$ data frame with named columns
#' @slot pathwaySet A list of known gene pathways
#' @slot response A numeric vector of length $N$: the dependent variable in a
#'   regression exercise
#'
#' @importFrom methods new
#'
#' @include createClass_OmicsPath.R
#' @seealso \code{"\link[=OmicsPathway-class]{OmicsPathway}"}
#'
#' @export

create_OmicsReg <- setClass("OmicsReg",
                            slots = c(response = "numeric"),
                            contains = "OmicsPathway")
