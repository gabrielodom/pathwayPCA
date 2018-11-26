#' Subset Pathway-Specific Data
#'
#' @description Given an \code{Omics} object and the name of a pathway, return
#'    the -omes in the assay and the response as a (tibble) data frame.
#'
#' @param object An object of class \code{OmicsPathway}, or an object extending
#'    this class.
#' @param pathName The name of a pathway contained in the pathway collection in
#'    the object.
#' @param ... Dots for additional internal arguments (currently unused).
#'
#' @return A data frame of the columns of the assay in the \code{Omics} object
#'    which are listed in the specified pathway. If the \code{Omics} object has
#'    response information, these are also included as the first column(s) of
#'    the data frame. If you have the suggested \code{tidyverse} package suite
#'    loaded, then this data frame will print as a \code{\link[tibble]{tibble}}.
#'    Otherwise, it will print as a data frame.
#'
#'
#' @details This function subsets the assay by the matching gene symbols or IDs
#'    in the specified pathway.
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#' @include createClass_OmicsSurv.R
#' @include createClass_OmicsReg.R
#' @include createClass_OmicsCateg.R
#' @include accessClass_OmicsPath.R
#' @include accessClass_OmicsRegCateg.R
#' @include accessClass_OmicsSurv.R
#'
#' @importFrom methods setGeneric
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   colon_Omics <- CreateOmics(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     response = colonSurv_df[, 1:2],
#'     respType = "survival"
#'   )
#'
#'   SubsetPathwayData(
#'     colon_Omics,
#'     "KEGG_RETINOL_METABOLISM"
#'   )
#'
#' }
#'
#' @rdname SubsetPathwayData
setGeneric("SubsetPathwayData",
           function(object, pathName, ...){
             standardGeneric("SubsetPathwayData")
           })


#' @rdname SubsetPathwayData
setMethod(f = "SubsetPathwayData", signature = "OmicsPathway",
          definition = function(object, pathName, ...){

            ###  The Design Matrix  ###
            path_ls <- getTrimPathwayCollection(object)[[pathName]]
            genes_char <- path_ls$IDs
            design_df <- getAssay(object)[, genes_char]

            ###  The Response Vectors  ###
            obj_class <- class(object)
            switch(obj_class,
                   OmicsPathway = {
                     out_df <- design_df
                   },
                   OmicsSurv = {
                     out_df <- data.frame(
                       EventTime = getEventTime(object),
                       EventObs  = getEvent(object),
                       design_df
                     )
                   },
                   OmicsReg = {
                     out_df <- data.frame(
                       Response = getResponse(object),
                       design_df
                     )
                   },
                   OmicsCateg = {
                     out_df <- data.frame(
                       Response = getResponse(object),
                       design_df
                     )
                   }
            )

            ###  Return  ###
            class(out_df) <- c("tbl_df", "tbl", "data.frame")
            out_df

          })
