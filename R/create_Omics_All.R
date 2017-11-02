#' Generation Functions for -Omics.* Class Objects
#'
#' These functions create valid objects of the "OmicsPathway", "OmicsSurv",
#'   "OmicsReg", and "OmicsCateg" classes.
#'
#' @section OmicsPathway:
#' Valid OmicsPathawy objects will have no response information, just the mass
#'   spectrometry (gene "design") matrix and the pathway list. OmicsPathway
#'   objects should be created only when unsupervised pathway extraction is
#'   needed. Because of the missing response, no pathway testing can be
#'   performed on an OmicsPathway object.
#'
#' @param massSpec_df An $N x p$ data frame with named columns
#' @param pathwaySet_ls A list of known gene pathways
#' @return An object of classes "OmicsPathway", "OmicsSurv", "OmicsReg", or
#'   "OmicsCateg".
#'
#' @include validClass_Omics.R
#' @include createClass_OmicsPath.R
#' @include createClass_OmicsSurv.R
#' @include createClass_OmicsReg.R
#' @include createClass_OmicsCateg.R
#'
#' @importFrom methods new
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsPath <- function(massSpec_df, pathwaySet_ls){

  new("OmicsPathway",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls)

}

#' @section OmicsSurv:
#' Valid OmicsSurv objects will have two response vectors: a vector of the most
#'   recently recorded follow-up times and a logical vector if that time marks
#'   an event (TRUE = observed event; FALSE = right-censored observation).
#'
#' @param eventTime_vec A numeric vector with $N$ observations corresponding to
#'   the last observed time of follow up
#' @param eventObserved_vec A logical vector with $N$ observations indicating
#'   right censoring. The values will be FALSE if the observation was censored
#'   (i.e., we did not observe an event).
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsSurv <- function(massSpec_df,
                             pathwaySet_ls,
                             eventTime_vec,
                             eventObserved_vec){

  new("OmicsSurv",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls,
      eventTime = eventTime_vec,
      eventObserved = eventObserved_vec)

}




#' @section OmicsReg and OmicsCateg:
#' Valid OmicsReg and OmicsCateg objects with have one response vector of
#'   continuous or categorial (as a factor) observations, respectively.
#'
#' @param response_num A numeric vector of length $N$: the dependent variable
#'   in a regression exercise
#'
#' @export
#' @rdname create_OmicsPathway
create_OmicsReg <- function(massSpec_df,
                            pathwaySet_ls,
                            response_num){

  new("OmicsReg",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls,
      response = response_num)

}

#' @param response_fact A factor vector of length $N$: the dependent variable of a
#'   generalized linear model
#' @export
#' @rdname create_OmicsPathway
create_OmicsCateg <- function(massSpec_df,
                              pathwaySet_ls,
                              response_fact){

  new("OmicsCateg",
      massSpec = massSpec_df,
      pathwaySet = pathwaySet_ls,
      response = response_fact)

}
