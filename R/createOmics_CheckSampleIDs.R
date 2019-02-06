#' Check Input Sample IDs
#'
#' @description Check the class of the sample IDs and if they are unique. This
#'    assumes that the sample IDs are in the first column.
#'
#' @param df An assay or phenotype data frame supplied to the
#'    \code{\link{CreateOmics}} function
#'
#' @return The same data frame, \strong{if} the sample IDs pass sanity checks,
#'    with the sample IDs as a character vector.
#'
#' @details This function checks that the sample IDs are unique, then coerces
#'    them from factor to character (if necessary), stores these IDs as the
#'    first column, then returns the same data frame.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'  # DO NOT CALL THIS FUNCTION DIRECTLY. CALL FROM WITHIN CreateOmics().
#'
#'  data("colonSurv_df")
#'  CheckSampleIDs(colonSurv_df[, -(2:3)])
#'
CheckSampleIDs <- function(df){
  # browser()

  dfSamp <- df[, 1, drop = TRUE]

  ###  Check for Row Names  ###
  dfRowNames <- row.names(df)
  genericRowNames <- as.character(seq_len(nrow(df)))
  defaultNames <- is.null(dfRowNames) || identical(dfRowNames, genericRowNames)
  if(!defaultNames){
    warning("Row names will be ignored. Sample IDs must be in the first column of the data
  frame. \n",
      immediate. = TRUE
    )
  }

  ###  Check Classes  ###
  if(inherits(dfSamp, "factor")){
    dfSamp_char <- levels(dfSamp)[dfSamp]
  } else if(inherits(dfSamp, "character")){
    dfSamp_char <- dfSamp
  } else {
    stop("Sample IDs must be stored in the first column of the assay or phenotype data
  frame. These IDs must be or extend class factor or character.")
  }

  ###  Check for Uniqueness / Missingness  ###
  if(anyNA(dfSamp_char)){
    stop("Missing sample IDs are not allowed.")
  }
  if(any(duplicated(dfSamp_char))){
    stop("Sample IDs must be unique.")
  }

  df[, 1] <- dfSamp_char

  ###  Return  ###
  df

}
