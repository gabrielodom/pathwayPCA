# Test CreateOmics wrapper
# source Test_ResponseConversion.R first

# CreateOmics <- function(assayData_df,
#                         sampleIDs_char,
#                         pathwayCollection_ls,
#                         response = NULL,
#                         respType = c("none", "survival", "regression", "categorical")){
#
#   # browser()
#
#   respType <- match.arg(respType)
#
#   if(respType != "none" && is.null(response)){
#     stop(paste0("Response must be specified for type = ", respType))
#   }
#   if(respType == "none" && !is.null(response)){
#     stop("Response type required when a response is given.")
#   }
#
#   if(!is.null(response)){
#     respClean <- convertResponse(response, type = respType)
#   }
#
#   switch (respType,
#     none = {
#
#       message("Creating object of class OmicsPathway.")
#       CreateOmicsPath(
#         assayData_df = assayData_df,
#         sampleIDs_char = sampleIDs_char,
#         pathwayCollection_ls = pathwayCollection_ls
#       )
#
#     },
#     survival = {
#
#       message("Creating object of class OmicsSurv.")
#       CreateOmicsSurv(
#         assayData_df = assayData_df,
#         sampleIDs_char = sampleIDs_char,
#         pathwayCollection_ls = pathwayCollection_ls,
#         eventTime_num = respClean$time,
#         eventObserved_lgl = respClean$dead
#       )
#
#     },
#     regression = {
#
#       message("Creating object of class OmicsReg.")
#       CreateOmicsReg(
#         assayData_df = assayData_df,
#         sampleIDs_char = sampleIDs_char,
#         pathwayCollection_ls = pathwayCollection_ls,
#         response_num = respClean
#       )
#
#     },
#     categorical = {
#
#       message("Creating object of class OmicsCateg.")
#       CreateOmicsCateg(
#         assayData_df = assayData_df,
#         sampleIDs_char = sampleIDs_char,
#         pathwayCollection_ls = pathwayCollection_ls,
#         response_fact = respClean
#       )
#
#     }
#   )
#
# }


library(pathwayPCA)

# # Test
# CreateOmics()
# CreateOmics(response = "12")
# CreateOmics(respType = "s")
# # Pathway
# CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
#              pathwayCollection_ls = gene_set_ls)
# # Survival
# CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
#              pathwayCollection_ls = gene_set_ls,
#              response = joinedExperiment_df[, 2:3],
#              respType = "s")
# # Regression
# CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
#              pathwayCollection_ls = gene_set_ls,
#              response = joinedExperiment_df[, 2],
#              respType = "r")
# # Categorical
# CreateOmics(assayData_df = joinedExperiment_df[, -c(1:3)],
#              pathwayCollection_ls = gene_set_ls,
#              response = joinedExperiment_df[, 3],
#              respType = "c")



gmt_path <- system.file(
  "extdata", "wikipathways_human_symbol.gmt",
  package = "pathwayPCA", mustWork = TRUE
)
wikipathways_PC <- read_gmt(gmt_path, description = TRUE)

data_path <- system.file(
  "extdata", "ovarian_PNNL_survival.RDS",
  package = "pathwayPCA", mustWork = TRUE
)
ovSurv_df <- readRDS(data_path)

ov_OmicsSurv <- CreateOmics(
  assayData_df = ovSurv_df[, -(2:3)],
  pathwayCollection_ls = wikipathways_PC,
  response = ovSurv_df[, 1:3],
  respType = "survival",
  minPathSize = 5,
  centerScale = c(TRUE, TRUE)
)

ov_OmicsSurv <- CreateOmics(
  assayData_df = ovSurv_df[, -(2:3)],
  pathwayCollection_ls = wikipathways_PC,
  response = ovSurv_df[-1, 1:3],
  respType = "survival",
  minPathSize = 5,
  centerScale = c(TRUE, TRUE)
)
