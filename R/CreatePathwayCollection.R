#' Manually Create a \code{pathwayCollection}-class Object.
#'
#' @description Manually create a \code{pathwayCollection} list similar to the
#'    output of the \code{\link{read_gmt}} function.
#'
#' @param sets_ls A named list of character vectors. Each vector should contain
#'    the names of the individual genes, proteins, sits, or CpGs within that
#'    set as a vector of character strings. If you create this pathway
#'    collection to integrate with data of class \code{Omics}*, the names
#'    contained in these vectors should have non-empty overlap with the feature
#'    names of the assay data frame that will be paired with this list in the
#'    subsequent analysis.
#' @param TERMS A character vector the same length as the \code{sets_ls} list
#'    with the proper names of the sets.
#' @param setType What is the type of the set: pathway set of gene, gene sites
#'    in RNA or DNA, or regions of CpGs. Defaults to \code{''pathway''}.
#' @param ... Additional vectors or data components related to the
#'    \code{sets_ls} list. These values should be passed as a name-value pair.
#'    See "Details" for more information.
#'
#' @return A list object with class \code{pathwayCollection}.
#'
#' @details This function checks the set list and set term inputs and then
#'    creates a \code{pathwayCollection} object from them. Pass additional
#'    list elements (such as the \code{description} of each set) using the
#'    form \code{tag = value} through the \code{...} argument (as in the
#'    \code{\link{list}} function). Because some functions in the
#'    \code{pathwayPCA} package add and edit elements of \code{pathwayCollection}
#'    objects, please do not create \code{pathwayCollection} list items named
#'    \code{setsize} or \code{n_tested}.
#'
#' @export
#'
#' @seealso \code{\link{read_gmt}}
#'
#' @examples
#'   data("colon_pathwayCollection")
#'
#'   CreatePathwayCollection(
#'     sets_ls = colon_pathwayCollection$pathways,
#'     TERMS = colon_pathwayCollection$TERMS
#'   )
#'
CreatePathwayCollection <- function(sets_ls, TERMS,
                                    setType = c("pathways", "genes", "regions"),
                                    ...){

  ###  Class Checks  ###
  if(!is.list(sets_ls)){
    stop("The sets_ls object must be a list of set membership vectors.")
  }

  if(!is.atomic(TERMS)){
    stop("The TERMS object must be an atomic vector.")
  }

  ###  Validation Checks  ###
  nPaths <- length(sets_ls)
  if(length(TERMS) != nPaths){
    stop("The TERMS vector must have the same length as the sets list.")
  }

  dotNames <- names(list(...))
  if(any(c("setsize", "n_tested") %in% dotNames)){
    warning("The names 'setsize' and 'n_tested' are reserved names of a pathwayCollection object.
            Values stored with these names may be overwritten during pathwayPCA function execution.
            Use with extreme caution.")
  }

  ###  Create and Return pathwayCollection  ###
  out_ls <- list(
    sets = sets_ls,
    TERMS = TERMS,
    ...
  )
  
  setType <- match.arg(setType)
  names(out_ls)[1] <- setType
  class(out_ls) <- c("pathwayCollection", "list")

  out_ls

}
