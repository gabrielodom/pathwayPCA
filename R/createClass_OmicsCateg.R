#' An S4 class for categorical responses within an \code{OmicsPathway} object
#'
#' @description This creates the \code{OmicsCateg} class which extends the
#'   \code{OmicsPathway} master class.
#'
#' @slot assayData_df An \eqn{N \times p} data frame with named columns.
#' @slot pathwaySet A list of known gene pathways with two elements:
#' \itemize{
#'   \item{\code{pathways} : }{A named list of character vectors. Each vector
#'      contains the names of the individual genes within that pathway as a
#'      vector of character strings. The names contained in these vectors must
#'      have non-empty overlap with the \emph{column names} of the
#'      \code{assayData_df} data frame. The names of the pathways (the list
#'      elements themselves) should be the a shorthand representation of the
#'      full pathway name.}
#'   \item{\code{TERMS} : }{A character vector the same length as the
#'     \code{pathways} list with the proper names of the pathways.}
#'   \item{\code{setsize} : }{A named integer vector the same length as the
#'     \code{pathways} list with the number of genes in each pathway. This list
#'     item is calculated during the creation step of a \code{create_OmicsCateg}
#'     function call.}
#' }
#' @slot response A \code{factor} vector of length \eqn{N}: the dependent
#'    variable of a generalized linear regression exercise. Currently, we
#'    support binary factors only. We expect to extend support to n-ary
#'    responses in the next package version.
#'
#' @importFrom methods setClass
#'
#' @include createClass_OmicsPath.R
#' @include createClass_validOmics.R
#'
#' @seealso \code{\link[=OmicsPathway-class]{OmicsPathway}},
#'   \code{\link{create_OmicsCateg}}
#'
#' @export
setClass("OmicsCateg",
         slots = c(response = "factor"),
         validity = valid_OmicsCateg,
         contains = "OmicsPathway")
