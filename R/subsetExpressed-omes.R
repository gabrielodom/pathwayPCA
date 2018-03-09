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
#' @param returnClass Should the returned object be of the same \code{Omics*}
#'   class as \code{object}, or should it be a list of pathway matrices? For
#'   internal calls from the \code{\link{extract_aesPCs}} function, this will be
#'   a \code{list} object of each pathway expressed gene matrix. For external
#'   calls, this defaults to the class of the object, allowing the user to input
#'   an object of class \code{OmicsSurv}, \code{OmicsReg}, \code{OmicsCateg},
#'   and have an object of the same class returned.
#' @param ... Dots for additional internal arguments (as necessary)
#'
#' @return If \code{returnClass = class(object)}: A valid \code{Omics*}-class
#'   object. This output object will be identical to the input object, except
#'   that any genes present in the pathways list, but not present in the MS
#'   design matrix, will have been removed. Additionally, the pathway list will
#'   have the number of genes in each trimmed pathway stored as the
#'   \code{trim_setsize} object.
#'
#'   If \code{returnClass = "list"}: A list of pathways with -ome measures
#'   expressed in the MS design matrix. Each element of the list will be named
#'   by its pathway, and the elements will be subset matrices of the original
#'   MS design matrix. See "details" for more information.
#'
#' @details This function takes in a data frame with named columns and a pathway
#'   list, all through one of the \code{Omics.*} classes. This function will
#'   then iterate over the list of pathways, extract columns from the MS design
#'   matrix which match the genes listed in that pathway, and remove any
#'   pathways with fewer than \code{trim} expressed genes. These matrices are
#'   returned as a named list (if \code{list} output is requested); or the genes
#'   not expressed in the MS design matrix are removed from the pathway list (if
#'   \code{class(object)} output is requested).
#'
#'   Note that some genes will be included in more than one pathway, so these
#'   pathways are not mutually exclusive. Further note that there may be many
#'   genes in the MS design matrix that are not included in the pathway sets, so
#'   these will not be extracted to the list. It is then vitally important to
#'   use either a very broad and generic pathway set list or a pathway set list
#'   that is appropriate for the mass spectrometry data supplied.
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @examples
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colonGenesets_ls")
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_OmicsSurv(massSpec_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colonGenesets_ls,
#'                                       eventTime_vec = colonSurv_df$OS_time,
#'                                       eventObserved_vec = as.logical(colonSurv_df$OS_event))
#'
#'   ###  Extract Expressed Genes  ###
#'   expressedOmes(colon_OmicsSurv)
#'
#'
#' @importFrom methods setGeneric
#' @rdname expressedOmes
setGeneric("expressedOmes",
           function(object,
                    trim = 3,
                    returnClass = class(object),
                    ...){
             standardGeneric("expressedOmes")
           }
)

#' @rdname expressedOmes
setMethod(f = "expressedOmes", signature = "OmicsPathway",
          definition = function(object,
                                trim = 3,
                                returnClass = class(object),
                                ...){

            # browser()

            genelist <- colnames(object@massSpec)
            paths_ls <- object@pathwaySet$pathways
            genesInPathway_vec <- unique(do.call(c, paths_ls))
            trimSetsize <- object@pathwaySet$setsize

            # Delete the genes from the pathways if they aren't recorded in our
            #   data matrix
            newPaths <- sapply(paths_ls, function(x){
              x[x %in% genelist]
            })

            trimSetsize <- sapply(seq_along(trimSetsize), function(i){

              size <- length(newPaths[[i]])
              ifelse(size < trim, 0, size)

            })
            names(trimSetsize) <- names(object@pathwaySet$setsize)

            # Remove any pathway that now has fewer than "trim" genes
            newPaths_trim <- sapply(seq_along(newPaths), function(i){

              shortPath <- length(newPaths[[i]]) < trim
              if(shortPath){
                NULL
              } else {
                newPaths[[i]]
              }

            })

            names(newPaths_trim) <- names(paths_ls)
            genesInTrimPathway_vec <- unique(do.call(c, newPaths_trim))
            # # If nullPaths resolves to int(0) (because the which() function
            # #   didn't find anything), then subsetting by -nullPaths will
            # #   give us an empty list. Filter(), however, can handle NULL and
            # #   integer(0) / character(0) list entries. We still need the
            # #   nullPaths object for later.
            nullPaths <- which(sapply(newPaths_trim, is.null))
            # paths_ls <- newPaths_trim[-nullPaths]
            cleanPaths_ls <- Filter(length, newPaths_trim)

            ###  Reporting Gene Overlap  ###
            missingPaths_char <- names(newPaths_trim)[nullPaths]
            pRmFeatures <- 1 -
              length(genesInTrimPathway_vec) / length(genesInPathway_vec)
            pSelectFeatures <- length(genesInTrimPathway_vec) / length(genelist)

            ###  Create the Return Object  ###
            if(as.character(returnClass) == "list"){

              extractedMatrix_ls <- lapply(cleanPaths_ls, function(x){
                object@massSpec[x]
              })


              attr(extractedMatrix_ls,
                   "missingPaths") <- missingPaths_char
              attr(extractedMatrix_ls,
                   "selectedFeature%") <- pSelectFeatures * 100
              attr(extractedMatrix_ls,
                   "removedFeature%") <- pRmFeatures * 100
              out <- extractedMatrix_ls

            } else {

              ###  Print Messages  ###
              # To display at 80 characters on the screen, our message can run
              #   from column 30 to 110 in te code.
              message(sprintf("Of the %i unique genes in the input pathway set, %.1f%% were not expressed in
  the input data and were therefore removed.",
                              length(genesInPathway_vec), pRmFeatures * 100))

              message(sprintf("After trimming unexpressed genes from the %i supplied pathways, we removed %i
  pathway(s) because they contained %i or fewer genes.",
                              length(paths_ls), length(missingPaths_char), trim))

              message(sprintf("Of the %i measured genes in the input data frame, %.1f%% were included in at
  least one pathway after trimming. \n",
                              length(genelist), pSelectFeatures * 100))

              attr(cleanPaths_ls,
                   "missingPaths") <- missingPaths_char
              attr(cleanPaths_ls,
                   "selectedFeature%") <- pSelectFeatures * 100
              attr(cleanPaths_ls,
                   "removedFeature%") <- pRmFeatures * 100
              out <- object
              out@pathwaySet$pathways <- cleanPaths_ls
              out@pathwaySet$trim_setsize <- trimSetsize

            }

            out

          })


# # Test:

# library(methods)
# library(pathwayPCA)
# load("data/ovarianFiltered_df.rda")
# load("data/genesets_ls.rda")
#
#
# testRedPath <- expressedOmes(testOmicsPath)
# testRedPath2 <- expressedOmes(testOmicsSurv)
# testRedPath3 <- expressedOmes(testOmicsReg)
# testRedPath4 <- expressedOmes(testOmicsCateg)
# identical(testRedPath, testRedPath2)
# identical(testRedPath, testRedPath3)
# identical(testRedPath, testRedPath4)
# # It works
