#' Rank the Top Genes from a Ranked Pathways Data Frame
#'
#' @description Given a supervised \code{Omics.*}-class object and a ranked
#'   pathways data frame returned by either the \code{\link{AESPCA_pVals}} or
#'   \code{\link{superPCA_pVals}} functions, rank the genes / proteins / lipids
#'   / metabolomes contained in each gene pathway by the weighted significance
#'   of their container pathways.
#'
#' @param object An object of class \code{OmicsSurv}, \code{OmicsReg}, or
#'   \code{OmicsCateg}.
#' @param pVals_df The ranked pathways data frame returned by either the
#'   \code{AESPCA_pVals} or \code{superPCA_pVals} functions.
#' @param percentile Return the most significant percent of the features
#'   contained in all pathways. Defaults to 0.05.
#'
#'
#' @return A named numeric vector. The names are the genes, and the values are
#'   the scores for those genes.
#'
#' @details This function takes in the pathway set information in a valid
#'   \code{Omics*}-class object and a data frame of ranked pathways (as returned
#'   by one of the two \code{*PCA_pVals()} functions). This function creates a
#'   matrix with pathways as the columns and all genes included in those
#'   pathways as the rows: the \eqn{i, j} entry of the matrix equals 1 if gene
#'   \eqn{i} is an element of pathway \eqn{j}. This is created after trimming
#'   the pathways to the assay data frame supplied using the
#'   \code{\link{expressedOmes}} function). The \code{topGenes} function then
#'   multiplies each pathway membership indicator column by the negative natural
#'   logarithm of the adjusted \eqn{p}-values for that pathway; if multiple FDR
#'   adjustment methods are used, then the score is the average of each negative
#'   logged \eqn{p}-value. This function then returns a named numeric vector of
#'   the sums of these scores for each gene, sorted in descending order.
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
#'   colon_OmicsSurv <- create_OmicsSurv(assayData_df = colonSurv_df[, -(1:2)],
#'                                       pathwaySet_ls = colonGenesets_ls,
#'                                       eventTime_vec = colonSurv_df$OS_time,
#'                                       eventObserved_vec = as.logical(colonSurv_df$OS_event))
#'
#'   ###  Extract PCs  ###
#'   surv_pVals_df <- superPCA_pVals(object = colon_OmicsSurv,
#'                                   parallel = TRUE,
#'                                   numCores = 2,
#'                                   adjustpValues = TRUE,
#'                                   adjustment = c("Hoch", "SidakSD"))
#'
#'   ###  Rank Genes  ###
#'   topGenes(object = colon_OmicsSurv,
#'            pVals_df = surv_pVals_df)
#'
#'
#' @importFrom methods setGeneric
#' @importFrom stats quantile
#' @rdname topGenes
setGeneric("topGenes",
           function(object, pVals_df, percentile = 0.05){
             standardGeneric("topGenes")
           }
)

#' @rdname topGenes
setMethod(f = "topGenes", signature = "OmicsPathway",
          definition = function(object, pVals_df, percentile = 0.05){

            # browser()

            clean_obj <- expressedOmes(object, message = FALSE)
            paths_ls <- clean_obj@pathwaySet$pathways
            genes_char <- unique(do.call(c, paths_ls))
            rm(object, clean_obj)


            ###  Calculate Pathway Scores  ###
            ranks_df <- pVals_df[, 1, drop = FALSE]
            # Currently, the first five columns are pathways, setsize, trim_size, terms,
            #   and rawp. The adjusted p-value columns do not start until column 6.
            adj_pVals_df <- pVals_df[, 6:ncol(pVals_df), drop = FALSE]
            # Add in a buffer in case one of the p-values is identically 0 (this is
            #   probable if the number of permutations is 1000 or smaller).
            minp <- min(unlist(adj_pVals_df)[unlist(adj_pVals_df) > 0])
            ranks_df$score <- rowMeans(-log(adj_pVals_df + minp / 100))

            ###  Create a Matrix of Pathway Membership  ###
            membership_mat <- sapply(1:nrow(pVals_df), function(path){

              idx <- which(names(paths_ls) == pVals_df[path, 1])
              (genes_char %in% paths_ls[[idx]]) * ranks_df[path, 2]

            })
            rownames(membership_mat) <- genes_char


            ###  Return the Top 5% Most Significant Genes  ###
            geneScores_num <- rowSums(membership_mat)
            cutoff <- quantile(geneScores_num, 1 - percentile)
            sort(geneScores_num[geneScores_num > cutoff], decreasing = TRUE)

          })

