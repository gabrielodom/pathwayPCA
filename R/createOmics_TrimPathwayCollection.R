#' Delete -Ome symbols or IDs without matching features recorded in a given
#'    assay data frame from a pathway collection
#'
#' @description Given a bio-assay design matrix and a \code{pathwayCollection}
#'    gene pathways list (each within an \code{Omics*}-class object), delete
#'    the genes / proteins / lipids / metabolomes / transcriptomes symbols or
#'    IDs recorded in each pathway which are not recorded in the assay data
#'    frame.
#'
#' @param object An object of class \code{OmicsPathway}, \code{OmicsSurv},
#'    \code{OmicsReg}, or \code{OmicsCateg}.
#' @param trim The minimum cutoff of matching -Ome measures before a pathway
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
#'    each trimmed pathway stored as the \code{n_tested} object.
#'
#' @details This function takes in a data frame with named columns and a
#'    \code{pathwayCollection} list, all through one of the \code{Omics*}
#'    classes. This function will then copy the pathway collection, iterate over
#'    the list of copied pathways, delete symbols or IDs from that pathway
#'    without matches from the bio-assay design matrix column names, and remove
#'    any pathways that have fewer than \code{trim} genes with corresponding
#'    columns in the assay. The genes not recorded in the bio-assay design
#'    matrix are removed from the copy of the pathway collection (the
#'    \code{trimPathwayCollection} object), but remain in the original pathway
#'    collection.
#'
#'    NOTE: some genes will be included in more than one pathway, so these
#'    pathways are not mutually exclusive. Further note that there may be many
#'    genes in the assay design matrix that are not included in the pathway
#'    sets, so these will not be extracted to the list. It is then vitally
#'    important to use either a very broad and generic \code{pathwayCollection}
#'    list or a \code{pathwayCollection} list that is appropriate for the assay
#'    data supplied. While you can create your own pathway lists, create proper
#'    \code{pathwayCollection} list objects by importing \code{.gmt} files with
#'    the \code{\link{read_gmt}} function.
#'
#' @keywords internal
#'
#' @export
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @examples
#' # DO NOT CALL THIS FUNCTION DIRECTLY. USE CreateOmics() INSTEAD.
#'
#'
#' @importFrom methods setGeneric
#' @rdname IntersectOmicsPwyCollct
setGeneric("IntersectOmicsPwyCollct",
           function(object,
                    trim = 3,
                    message = TRUE,
                    ...){
             standardGeneric("IntersectOmicsPwyCollct")
           }
)

#' @rdname IntersectOmicsPwyCollct
setMethod(f = "IntersectOmicsPwyCollct", signature = "OmicsPathway",
          definition = function(object,
                                trim = 3,
                                message = TRUE,
                                ...){

            # browser()

            genelist <- colnames(object@assayData_df)
            paths_ls <- object@pathwayCollection$pathways
            genesInPathway_vec <- unique(do.call(c, paths_ls))

            # Delete the genes from the pathways if they aren't recorded in our
            #   data matrix
            newPaths <- sapply(paths_ls, function(x){
              x[x %in% genelist]
            })

            newPathsLen <- lengths(newPaths)
            trimSetsize <- ifelse(newPathsLen < trim, 0, newPathsLen)

            # Remove any pathway that now has fewer than "trim" genes
            newPaths_trim <- sapply(seq_along(newPaths), function(i){

              if(trimSetsize[i] == 0){
                NULL
              } else {
                newPaths[[i]]
              }

            })

            names(newPaths_trim) <- names(paths_ls)
            cleanPaths_ls <- Filter(length, newPaths_trim)

            # If nullPaths resolves to int(0) (because the which() function
            #   didn't find anything), then subsetting by -nullPaths will give
            #   us an empty list. Filter(), however, can handle NULL and
            #   integer(0) / character(0) list entries, which is why we use it
            #   above.
            nullPaths <- which(sapply(newPaths_trim, is.null))
            TERMS_char <- object@pathwayCollection$TERMS
            description_char <- object@pathwayCollection$description
            setsize_int <- object@pathwayCollection$setsize
            if(length(nullPaths) > 0){

              newTERMS_char <- TERMS_char[-nullPaths]

              if(is.null(description_char)){
                newDesc_char <- description_char
              } else {
                newDesc_char <- description_char[-nullPaths]
              }

              newSetsize_int <- setsize_int[-nullPaths]
              trimSetsize <- trimSetsize[-nullPaths]

            } else {

              newTERMS_char <- TERMS_char
              newDesc_char <- description_char
              newSetsize_int <- setsize_int

            }

            missingPaths_char <- names(newPaths_trim)[nullPaths]


            ###  Reporting Gene Overlap  ###
            genesInTrimPathway_vec <- unique(do.call(c, newPaths_trim))
            pRmFeatures <- 1 -
              length(genesInTrimPathway_vec) / length(genesInPathway_vec)
            pSelectFeatures <- length(genesInTrimPathway_vec) / length(genelist)

            ###  Print Messages  ###
            # To display at 80 characters on the screen, our message can run
            #   from column 30 to 110 in the code.
            if(message){

              message(
                sprintf("Of the %i unique genes in the input pathways list, %.1f%% were not found in the
  input assay data and were therefore removed.",
                  length(genesInPathway_vec), pRmFeatures * 100
                )
              )

              message(
                sprintf("After trimming these genes from the %i supplied pathways, we removed %i
  pathway(s) because they contained %i or fewer genes.",
                  length(paths_ls), length(missingPaths_char), trim
                )
              )

              message(
                sprintf("Of the %i measured genes in the input assay data, %.1f%% (%i) were included in at
  least one pathway after trimming. \n",
                  length(genelist), pSelectFeatures * 100,
                  length(genesInTrimPathway_vec)
                )
              )

            }

            if(length(genesInTrimPathway_vec) < 1){
              stop(
  "The feature names of the assay do not match the supplied features within the
pathway collection. Please supply a different pathway collection."
              )
            }


            ###  Create the Trimmed Pathway Collection  ###
            attr(cleanPaths_ls, "minFeatures") <- trim
            attr(cleanPaths_ls, "missingPaths") <- missingPaths_char
            attr(cleanPaths_ls, "selectedFeature%") <- pSelectFeatures * 100
            attr(cleanPaths_ls, "removedFeature%") <- pRmFeatures * 100

            trim_PC <- object@pathwayCollection
            trim_PC$pathways <- cleanPaths_ls
            trim_PC$TERMS <- newTERMS_char
            trim_PC$description <- newDesc_char
            trim_PC$setsize <- newSetsize_int
            trim_PC$n_tested <- trimSetsize


            ###  Create the Return Object  ###
            out <- object
            out@trimPathwayCollection <- trim_PC

            out

          })

