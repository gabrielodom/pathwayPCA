#' Extract expressed -Omes matching a pathways list from an assay data frame
#'
#' @description Given a bio-assay design matrix and a \code{pathwaySet} gene
#'    pathways list (each within an \code{Omics*}-class object), extract the
#'    genes / proteins / lipids / metabolomes / transcriptomes contained in each
#'    pathways list which are expressed in the assay data frame.
#'
#' @param object An object of class \code{OmicsPathway}, \code{OmicsSurv},
#'    \code{OmicsReg}, or \code{OmicsCateg}.
#' @param trim The minimum cutoff of expressed -Ome measures before a pathway
#'    is excluded. Defaults to 3.
#' @param message Should this function return diagnostic messages? Messages
#'    concern the percentage of genes included in the pathways list but not
#'    measured in the data, genes measured in the data but not called for in the
#'    pathways, and the number of pathways ignored due to too few number of
#'    genes present after trimming. Defaults to \code{TRUE}.
#' @param ... Dots for additional internal arguments (as necessary).
#'
#' @return A valid \code{Omics*}-class object. This output object will be
#'    identical to the input object, except that any genes present in the
#'    pathways list, but not present in the MS design matrix, will have been
#'    removed. Additionally, the pathway list will have the number of genes in
#'    each trimmed pathway stored as the \code{trim_setsize} object.
#'
#' @details This function takes in a data frame with named columns and a
#'    \code{pathwaySet} list, all through one of the \code{Omics*} classes.
#'    This function will then iterate over the list of pathways, extract columns
#'    from the bio-assay design matrix which match the genes listed in that
#'    pathway, and remove any pathways with fewer than \code{trim} expressed
#'    genes. The genes not expressed in the bio-assay design matrix are removed
#'    from the \code{pathwaySet} list.
#'
#'    NOTE: some genes will be included in more than one pathway, so these
#'    pathways are not mutually exclusive. Further note that there may be many
#'    genes in the assay design matrix that are not included in the pathway
#'    sets, so these will not be extracted to the list. It is then vitally
#'    important to use either a very broad and generic \code{pathwaySet} list
#'    or a \code{pathwaySet} list that is appropriate for the assay data
#'    supplied. While you can create your own pathway lists, create proper
#'    \code{pathwaySet} list objects by importing \code{.gmt} files with the
#'    \code{\link{read_gmt}} function.
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @examples
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwaySet")
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- create_OmicsSurv(assayData_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colon_pathwaySet,
#'                                       eventTime_num = colonSurv_df$OS_time,
#'                                       eventObserved_lgl = as.logical(colonSurv_df$OS_event))
#'
#'   ###  Extract Expressed Genes  ###
#'   expressedOmes(colon_OmicsSurv)
#' }
#'
#' @importFrom methods setGeneric
#' @rdname expressedOmes
setGeneric("expressedOmes",
           function(object,
                    trim = 3,
                    message = TRUE,
                    ...){
             standardGeneric("expressedOmes")
           }
)

#' @rdname expressedOmes
setMethod(f = "expressedOmes", signature = "OmicsPathway",
          definition = function(object,
                                trim = 3,
                                message = TRUE,
                                ...){

            # browser()

            genelist <- colnames(object@assayData_df)
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

            ###  Print Messages  ###
            # To display at 80 characters on the screen, our message can run
            #   from column 30 to 110 in the code.
            if(message){

              message(sprintf("Of the %i unique genes in the input pathways list, %.1f%% were not expressed in
  the input data and were therefore removed.",
                              length(genesInPathway_vec), pRmFeatures * 100))

              message(sprintf("After trimming unexpressed genes from the %i supplied pathways, we removed %i
  pathway(s) because they contained %i or fewer genes.",
                              length(paths_ls), length(missingPaths_char), trim))

              message(sprintf("Of the %i measured genes in the input data frame, %.1f%% were included in at
  least one pathway after trimming. \n",
                              length(genelist), pSelectFeatures * 100))

            }


            ###  Create the Return Object  ###
            # if(as.character(returnClass) == "list"){
            #
            #   extractedMatrix_ls <- lapply(cleanPaths_ls, function(x){
            #     object@assayData_df[x]
            #   })
            #
            #
            #   attr(extractedMatrix_ls,
            #        "missingPaths") <- missingPaths_char
            #   attr(extractedMatrix_ls,
            #        "selectedFeature%") <- pSelectFeatures * 100
            #   attr(extractedMatrix_ls,
            #        "removedFeature%") <- pRmFeatures * 100
            #   out <- extractedMatrix_ls
            #
            # } else {

              attr(cleanPaths_ls,
                   "missingPaths") <- missingPaths_char
              attr(cleanPaths_ls,
                   "selectedFeature%") <- pSelectFeatures * 100
              attr(cleanPaths_ls,
                   "removedFeature%") <- pRmFeatures * 100
              out <- object
              out@pathwaySet$pathways <- cleanPaths_ls
              out@pathwaySet$trim_setsize <- trimSetsize

            # }

            out

          })

