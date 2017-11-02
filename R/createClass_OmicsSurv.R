#' An S4 class for survival responses within an "OmicsPathway" object
#'
#' @description This creates the "OmicsSurv" class which extends the
#'   "OmicsPathway" master class.
#'
#' @slot massSpec An $N x p$ data frame with named columns
#' @slot pathwaySet A list of known gene pathways
#' @slot eventTime A numeric vector with $N$ observations corresponding to the
#'   last observed time of follow up
#' @slot eventObserved A logical vector with $N$ observations indicating right
#'   censoring. The values will be FALSE if the observation was censored (i.e.,
#'   we did not observe an event).
#'
#' @importFrom methods new
#'
#' @include createClass_OmicsPath.R validClass_Omics.R
#' @seealso \code{"\link[=OmicsPathway-class]{OmicsPathway}"}
#'
#' @export

create_OmicsSurv <- setClass("OmicsSurv",
                             slots = c(eventTime = "numeric",
                                       eventObserved = "logical"),
                             validity = valid_OmicsSurv,
                             contains = "OmicsPathway")
