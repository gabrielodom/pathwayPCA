#' Extract expressed -omes matching a gene set from a mass spectrometry matrix
#'
#' @description Given a mass spectrometry design matrix and a gene pathways
#'   list (each within an \code{Omics.*}-class object), extract the genes /
#'   proteins / lipids / metabolomes contained in each gene pathway set which
#'   are expressed in the MS design matrix.
#'
#' @param object An object of class \code{OmicsPathway}, \code{OmicsSurv},
#'   \code{OmicsReg}, or \code{OmicsCateg}
#' @param trim The minimum cutoff of expressed -ome measures before a pathway
#'   is excluded. Defaults to 3.
#' @param ... Dots for additional internal arguments (as necessary)
#'
#' @return A list of pathways with -ome measures expressed in the MS design
#'   matrix. Each element of the list will be named by its pathway, and the
#'   elements will be subset matrices of the original MS design matrix. See
#'   "details" for more information
#'
#' @details This function takes in a data frame with named columns and a pathway
#'   list, all through one of the \code{Omics.*} classes. This function will
#'   then iterate over the list of pathways, extract columns from the MS design
#'   matrix which match the genes listed in that pathway, and remove any
#'   pathways with fewer than \code{trim} expressed genes. These matrices are
#'   returned as a named list. Note thatsome genes will be included in more than
#'   one pathway, so these pathways are not mutually exclusive. Further note
#'   that there may be many genes in the MS design matrix that are not included
#'   in the pathway sets, so these will not be extracted to the list. It is then
#'   vitally important to use either a very broad and generic pathway set list
#'   or a pathway set list that is appropriate for the mass spectrometry data
#'   supplied.
#'
#' @export
#'
#' @include validClass_Omics.R
#' @include createClass_OmicsPath.R
#'
#' @importFrom methods setGeneric
#' @rdname expressedOmes
setGeneric("expressedOmes",
           function(object, trim = 3, ...){
             standardGeneric("expressedOmes")
           }
)

#' @rdname expressedOmes
setMethod(f = "expressedOmes", signature = "OmicsPathway",
          definition = function(object, trim = 3, ...){
            genelist <- colnames(object@massSpec)

            # Delete the genes from the pathways if they aren't recorded in our
            #   data matrix
            newPaths <- sapply(object@pathwaySet, function(x){
              x[x %in% genelist]
            })

            # Remove any pathway that now has fewer than "trim" genes
            newPaths_trim <- sapply(seq_along(newPaths), function(i){

              shortPath <- length(newPaths[[i]]) < trim
              if(shortPath){
                NULL
              } else {
                newPaths[[i]]
              }

            })

            names(newPaths_trim) <- names(object@pathwaySet)

            nullPaths <- which(sapply(newPaths_trim, is.null))
            paths_ls <- newPaths_trim[-nullPaths]

            extractedMatrix_ls <- lapply(paths_ls, function(x){
              object@massSpec[x]
            })


            attr(extractedMatrix_ls,
                 "missingPaths") <- names(newPaths_trim)[nullPaths]
            extractedMatrix_ls

          })
