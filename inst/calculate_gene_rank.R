#' Rank -Omes by adjusted significance given in a ranked-pathways data frame
#'
#' @description Given a supervised \code{Omics*}-class object and a ranked
#'    pathways data frame returned by either the \code{\link{AESPCA_pVals}} or
#'    \code{\link{SuperPCA_pVals}} functions, calculate the weighted rank the
#'    genes / proteins / lipids / metabolomes / transcriptomes contained in each
#'    pathway by the significance of their container pathways.
#'
#' @param object An object of class \code{OmicsSurv}, \code{OmicsReg}, or
#'    \code{OmicsCateg}.
#' @param pVals_df The ranked pathways data frame returned by either the
#'    \code{\link{AESPCA_pVals}} or \code{\link{SuperPCA_pVals}} functions.
#'    Missing \eqn{p}-values (from trimmed pathways) are omitted.
#' @param percentile Return the most significant \eqn{q} percent of the features
#'    contained in all pathways. Defaults to 0.01.
#'
#'
#' @return A list of two named numeric vectors. For both vectors, the names are
#'    the genes, and the values are the scores for those genes. The first vector
#'    is the sum of scores across all pathways; the second vector is this score
#'    sum divided by the number of pathways which contain that particular gene.
#'    The summed vector does not adjust for genes which appear more frequently
#'    in pathways, while the averaged vector does. See "Details" for more
#'    information.
#'
#' @details This function takes in the pathway set information in a valid
#'    \code{Omics*}-class object and a data frame of ranked pathways (as
#'    returned by either the \code{\link{AESPCA_pVals}} or
#'    \code{\link{SuperPCA_pVals}} functions). This function creates a matrix
#'    with pathways as the columns and all genes included in those pathways as
#'    the rows: the \eqn{i, j} entry of the matrix equals 1 if gene \eqn{i} is
#'    an element of pathway \eqn{j}. (This is created after trimming the
#'    pathways to the assay data frame supplied using the
#'    \code{\link{IntersectOmicsPwyCollct}} function.) The \code{topGenes}
#'    function then multiplies each pathway membership indicator column by the
#'    negative natural logarithm of the adjusted \eqn{p}-values for that
#'    pathway; if multiple FDR adjustment methods are used, then the score is
#'    the average of each negative logged \eqn{p}-values. This function then
#'    returns two named numeric vectors: the sum of these gene scores and the
#'    means of the non-zero gene scores, sorted in descending order.
#'
#' @include createClass_validOmics.R
#' @include createClass_OmicsPath.R
#'
#' @importFrom methods setGeneric
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#'   ###  Load the Example Data  ###
#'   data("colonSurv_df")
#'   data("colon_pathwayCollection")
#'
#'   ###  Create an OmicsSurv Object  ###
#'   colon_OmicsSurv <- CreateOmicsSurv(
#'     assayData_df = colonSurv_df[, -(1:2)],
#'     pathwayCollection_ls = colon_pathwayCollection,
#'     eventTime_num = colonSurv_df$OS_time,
#'     eventObserved_lgl = as.logical(colonSurv_df$OS_event)
#'   )
#'
#'   ###  Calculate Pathway p-Values  ###
#'   colonSurv_pVals_df <- SuperPCA_pVals(
#'     object = colon_OmicsSurv,
#'     parallel = TRUE,
#'     numCores = 2,
#'     adjustpValues = TRUE,
#'     adjustment = c("Hoch", "SidakSD")
#'   )
#'
#'   ###  Find the Top Genes  ###
#'   topGenes(object = colon_OmicsSurv, pVals_df = colonSurv_pVals_df)
#' }
#'
#' @rdname topGenes
setGeneric("topGenes",
           function(object, pVals_df, percentile = 0.01){
             standardGeneric("topGenes")
           }
)

#' @rdname topGenes
setMethod(f = "topGenes", signature = "OmicsPathway",
          definition = function(object, pVals_df, percentile = 0.01){

            # browser()

            clean_obj <- IntersectOmicsPwyCollct(object, message = FALSE)
            paths_ls <- clean_obj@pathwayCollection$pathways
            genes_char <- unique(do.call(c, paths_ls))
            rm(object, clean_obj)

            # Because pVals_df is probably a tibble, df[1, 1] will not print
            #   unless the tidyverse is loaded. Therefore, we need to force
            #   pVals_df to be a data frame. Also, because trimmed paths will
            #   have NA p-values, we'll remove these from any ranking.
            pVals_df <- na.omit(as.data.frame(pVals_df))


            ###  Calculate Pathway Scores  ###
            ranks_df <- pVals_df[, 1, drop = FALSE]
            # Currently, the first five columns are pathways, setsize, trim_size,
            #   terms, and rawp. The adjusted p-value columns do not start until
            #   column 6.
            adj_pVals_df <- pVals_df[, 6:ncol(pVals_df), drop = FALSE]
            # Add in a buffer in case one of the p-values is identically 0 (this
            #   is probable if the number of permutations is 1000 or smaller).
            adjP_vec <- unlist(adj_pVals_df)
            if(any(adjP_vec == 0)){
              minp <- min(adjP_vec[adjP_vec > 0])
            } else {
              minp <- 0
            }

            # By padding the 0 values, we can't take the log of zero. By padding
            #   the -log(1) values, we can't have a row sum of zero
            ranks_df$score <- -log(rowMeans(adj_pVals_df + minp / 100)) + 1e-99

            ###  Create a Matrix of Pathway Membership  ###
            membership_mat <- sapply(1:nrow(pVals_df), function(path){

              idx <- which(names(paths_ls) == pVals_df[path, 1])
              (genes_char %in% paths_ls[[idx]]) * ranks_df[path, 2]

            })
            rownames(membership_mat) <- genes_char


            ###  Return the Top 5% Most Significant Genes  ###
            # We will take the average non-zero score. If we had taken the sum,
            #   then a gene that showed up in a ton of unrelated pathways could
            #   have a higher score than a gene that shows up only once but in
            #   a strongly-related pathway. We can count the non-zero scores
            #   by using !!: the first turns all numbers into FALSE if the
            #   number does not equal 0 and TRUE otherwise, so we negate that
            #   again to turn all non-zero numbers into TRUE (which = 1).
            rawScore <- rowSums(membership_mat)
            rawCutoff <- quantile(rawScore, 1 - percentile)
            rawRank <- sort(rawScore[rawScore >= rawCutoff], decreasing = TRUE)

            adjScore <- rawScore / rowSums(!!membership_mat)
            adjCutoff <- quantile(adjScore, 1 - percentile)
            adjRank <- sort(adjScore[adjScore >= adjCutoff], decreasing = TRUE)

            list(summedRank = rawRank,
                 averagedRank = adjRank)

          })

