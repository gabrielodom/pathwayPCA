#'  Alternative Model t-Statistic Matrix
#'
#' @description A matrix of t-scores from a subset of the pathway models for
#'   the \code{supervised_Tumors_df} and \code{supervised_patInfo_df} data set.
#'   These pathways are a subset of the full pathway set where pathways smaller
#'   than 5 genes have been deleted.
#'
#' @format A 4,240 x 20 matrix of t-scores. Each row corresponds to the t-scores
#'   calculated when the patient survival information is fit to the 4,240 gene
#'   pathway sets. There are 20 t-scores, corresponding to 20 threshold values
#'   for the Cox Propotional Hazards regression model.
#'
#' @source Personal computation
"tScores4240_mat"
