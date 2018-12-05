#' Merge Phenotype and Assay Data by First Column (Sample ID)
#'
#' @description Match the records from the phenotype data to the values in the
#'    assay data by sample ID. Return rows from each data frame with matches in
#'    both data frames. The sample ID must be the first column in both data
#'    frames.
#'
#' @param pheno_df Phenotype data frame with the sample IDs in the first column
#' @param assay_df Assay data frame with the sample IDs in the first column
#'
#' @return A list of three elements:
#'   \itemize{
#'     \item{\code{assay} : }{A data frame with the rows from \code{assay_df}
#'       which are contained in \code{pheno_df}, ordered by their position in
#'       \code{pheno_df}.}
#'     \item{\code{response} : }{A data frame with the rows from
#'       \code{pheno_df} which are contained in \code{assay_df}.}
#'     \item{\code{sampleID} : }{A vector of the sample IDs shared by both data
#'       frames, ordered by their position in \code{pheno_df}.}
#'   }
#'
#' @details Don't use this function. This is simply a wrapper around the
#'    \code{\link{merge}} function with extra checks for the class of the ID
#'    column. If you want to merge your two data frames by sample ID, you should
#'    use the \code{inner_join} function from the \code{dplyr} package instead.
#'    It's easier. See \url{https://dplyr.tidyverse.org/reference/join.html}.
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'   # DO NOT CALL THIS FUNCTIONS DIRECTLY. USE CreateOmics() INSTEAD.
JoinPhenoAssay <- function(pheno_df, assay_df){

  ###  Check Phenotype Sample IDs  ###
  phSamp <- pheno_df[, 1, drop = TRUE]

  if(inherits(phSamp, "factor")){
    phSamp_char <- levels(phSamp)[phSamp]
  } else if(inherits(phSamp, "numeric")){
    phSamp_char <- as.character(phSamp)
  } else if(inherits(phSamp, "character")){
    phSamp_char <- phSamp
  } else {
    stop("Sample IDs should be stored in the first column of the phenotype data frame.
         These IDs must be or extend class numeric, factor, or character.")
  }

  pheno_df[, 1] <- phSamp_char

  ###  Check Assay Sample IDs
  assSamp <- assay_df[, 1, drop = TRUE]

  if(inherits(assSamp, "factor")){
    assSamp_char <- levels(assSamp)[assSamp]
  } else if(inherits(assSamp, "numeric")){
    assSamp_char <- as.character(assSamp)
  } else if(inherits(assSamp, "character")){
    assSamp_char <- assSamp
  } else {
    stop("Sample IDs should be stored in the first column of the assay data frame. These
         IDs must be or extend class numeric, factor, or character.")
  }

  assay_df[, 1] <- assSamp_char

  ###  Check Equality  ###
  keepIDs <- intersect(assSamp_char, phSamp_char)

  if(identical(assSamp_char, phSamp_char)){

    out_ls <- list(
      assay = assay_df[, -1],
      response = pheno_df[, -1, drop = FALSE],
      sampleID = phSamp_char
    )

  } else if(length(keepIDs) > 0){

    message(
      sprintf("There are %i samples shared by the assay and phenotype data.",
              length(keepIDs))
    )

    out_df <- merge(pheno_df, assay_df, 1)
    outClass_char <- union(
      class(pheno_df), class(assay_df)
    )
    class(out_df) <- outClass_char

    out_ls <- list(
      assay    = out_df[, -(1:ncol(pheno_df))],
      response = out_df[, 2:ncol(pheno_df), drop = FALSE],
      sampleID = out_df[, 1, drop = TRUE]
    )

  } else {
    stop("There are no samples with the same sample IDs in the assay and phenotype data.")
  }

  ###  Return  ###
  out_ls

  }
