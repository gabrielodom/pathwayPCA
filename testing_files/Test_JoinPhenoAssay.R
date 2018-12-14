# Test Joining Pheno and Assay Data
# Gabriel Odom
# 20181204

# We want to add the functionality to join the phenotype and assay data within
#   the CreateOmics() call. I strongly disagree with this, but I was overruled.
#   See Issue 42 (https://github.com/gabrielodom/pathwayPCA/issues/42) for more
#   information.

library(pathwayPCA)


######  Test Setup  ###########################################################
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



######  Joining Function  #####################################################
# JoinPhenoAssay <- function(pheno_df, assay_df){
#   browser()
#
#   ###  Check Phenotype Sample IDs  ###
#   phSamp <- pheno_df[, 1, drop = TRUE]
#
#   if(inherits(phSamp, "factor")){
#     phSamp_char <- levels(phSamp)[phSamp]
#   } else if(inherits(phSamp, "numeric")){
#     phSamp_char <- as.character(phSamp)
#   } else if(inherits(phSamp, "character")){
#     phSamp_char <- phSamp
#   } else {
#     stop("Sample IDs should be stored in the first column of the phenotype data frame.
#   These IDs must be or extend class numeric, factor, or character.")
#   }
#
#   pheno_df[, 1] <- phSamp_char
#
#   ###  Check Assay Sample IDs
#   assSamp <- assay_df[, 1, drop = TRUE]
#
#   if(inherits(assSamp, "factor")){
#     assSamp_char <- levels(assSamp)[assSamp]
#   } else if(inherits(assSamp, "numeric")){
#     assSamp_char <- as.character(assSamp)
#   } else if(inherits(assSamp, "character")){
#     assSamp_char <- assSamp
#   } else {
#     stop("Sample IDs should be stored in the first column of the assay data frame. These
#   IDs must be or extend class numeric, factor, or character.")
#   }
#
#   assay_df[, 1] <- assSamp_char
#
#   ###  Check Equality  ###
#   keepIDs <- intersect(assSamp_char, phSamp_char)
#
#   if(identical(assSamp_char, phSamp_char)){
#
#     out_ls <- list(
#       assay = assay_df[, -1],
#       response = pheno_df[, -1, drop = FALSE],
#       sampleID = phSamp_char
#     )
#
#   } else if(length(keepIDs) > 0){
#
#     message(
#       sprintf("There are %i samples shared by the assay and phenotype data.",
#               length(keepIDs))
#     )
#
#     out_df <- merge(pheno_df, assay_df, 1)
#     outClass_char <- union(
#       class(pheno_df), class(assay_df)
#     )
#     class(out_df) <- outClass_char
#
#     out_ls <- list(
#       assay    = out_df[, -(1:ncol(pheno_df))],
#       response = out_df[, 2:ncol(pheno_df), drop = FALSE],
#       sampleID = out_df[, 1]
#     )
#
#   } else {
#     stop("There are no samples with the same sample IDs in the assay and phenotype data.")
#   }
#
#   ###  Return  ###
#   out_ls
#
# }

# Test
JoinPhenoAssay(
  pheno_df = ovSurv_df[, 1:3],
  assay_df = ovSurv_df[, -(2:3)]
)

JoinPhenoAssay(
  pheno_df = ovSurv_df[-1, 1:3],
  assay_df = ovSurv_df[, -(2:3)]
)

JoinPhenoAssay(
  pheno_df = data.frame(
    IDs = sample(LETTERS, 83, replace = TRUE),
    ovSurv_df[, 2:3]
  ),
  assay_df = ovSurv_df[, -(2:3)]
)


# df1 <- data.frame(
#   IDs1 = as.character(1:6),
#   values = LETTERS[1:6],
#   stringsAsFactors = FALSE
# )
# df2 <- data.frame(
#   IDs2 = as.character(9:4),
#   values = 4:9,
#   stringsAsFactors = FALSE
# )
# intersect(df1$IDs1, df2$IDs2)
# intersect(df2$IDs2, df1$IDs1)
#
# # df1$IDs
# # df1$IDs %in% intersect(df1$IDs, df2$IDs)
# # df2$IDs
# # df2$IDs %in% intersect(df1$IDs, df2$IDs)
#
# merge(df1, df2, 1)
