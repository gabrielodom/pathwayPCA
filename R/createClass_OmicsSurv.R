#' An S4 class for survival responses within an \code{OmicsPathway} object
#'
#' @description This creates the \code{OmicsSurv} class which extends the
#'   \code{OmicsPathway} master class.
#'
#' @slot assayData_df An \eqn{N \times p} data frame with named columns.
#' @slot pathwayCollection A list of known gene pathways with three or four
#'    elements:
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
#'   \item{\code{description} : }{An optional character vector the same length
#'      as the \code{pathways} list with additional information about the
#'      pathways.}
#'   \item{\code{setsize} : }{A named integer vector the same length as the
#'     \code{pathways} list with the number of genes in each pathway. This list
#'     item is calculated during the creation step of a \code{CreateOmicsSurv}
#'     function call.}
#' }
#' @slot eventTime A \code{numeric} vector with \eqn{N} observations
#'   corresponding to the last observed time of follow up.
#' @slot eventObserved A \code{logical} vector with \eqn{N} observations
#'   indicating right-censoring. The values will be \code{FALSE} if the
#'   observation was censored (i.e., we did not observe an event).
#'
#' @importFrom methods setClass
#'
#' @include createClass_OmicsPath.R
#' @include createClass_validOmics.R
#'
#' @seealso \code{\link[=OmicsPathway-class]{OmicsPathway}},
#'   \code{\link{CreateOmicsSurv}}
#'
#' @export
setClass("OmicsSurv",
         slots = c(eventTime = "numeric",
                   eventObserved = "logical"),
         validity = ValidOmicsSurv,
         contains = "OmicsPathway")
