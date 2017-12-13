#'  Alternative Model t-Statistic Matrix
#'
#' @description A matrix of t-scores from pathway models for the
#'   \code{supervised_Tumors_df} and \code{supervised_patInfo_df} data set.
#'
#' @format A 7,949 x 20 matrix of t-scores. Each row corresponds to the t-scores
#'   calculated when the patient survival information is fit to the 7,949 gene
#'   pathway sets. There are 20 t-scores, corresponding to 20 threshold values
#'   for the Cox Propotional Hazards regression model.
#'
#' @source Personal computation
"tScores_mat"
