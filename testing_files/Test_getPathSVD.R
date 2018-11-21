# Draft of subsetting functions for PCs and Loadings
# Gabriel Odom
# 20181121


######  Setup  ################################################################
library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")

colon_OmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[, 1:2],
  respType = "survival"
)

colon_superpca <- SuperPCA_pVals(
  object = colon_OmicsSurv,
  numPCs = 1,
  parallel = TRUE,
  numCores = 16,
  adjustpValues = TRUE,
  adjustment = c("BH", "BY")
)

colon_aespca <- AESPCA_pVals(
  object = colon_OmicsSurv,
  numReps = 100,
  numPCs = 1,
  parallel = TRUE,
  numCores = 16,
  adjustpValues = TRUE,
  adjustment = c("BH", "BY")
)


######  getPathDecomp  ########################################################

getPathSVD <- function(pcOut, pathway_char, ...) UseMethod("getPathSVD")

getPathSVD.superpcOut <- function(pcOut, pathway_char){

  ###  Check for Matches  ###
  pathID_idx <- which(pcOut$pVals_df$pathways == pathway_char)
  term_idx   <- which(pcOut$pVals_df$terms == pathway_char)
  if(length(pathID_idx) + length(term_idx) == 0){
    stop("Supplied pathway does not match any pathway in the supplied object.")
  }

  ###  Match pathway ID to pathway name  ###
  if(length(term_idx) == 1){

    pathID  <- as.character(pcOut$pVals_df[term_idx, "pathways"])
    pathway <- pathway_char

  } else if(length(pathID_idx) == 1) {

    pathID  <- pathway_char
    pathway <- as.character(pcOut$pVals_df[pathID_idx, "terms"])

  }

  ###  Select the PCs and Loadings that Match this Pathway ID  ###
  list(
    PCs = pcOut$PCs_ls[[pathID]],
    Loadings = pcOut$loadings_ls[[pathID]],
    pathway = pathID,
    term = pathway
  )

}

getPathSVD.aespcOut <- getPathSVD.superpcOut


# Test
getPathSVD(colon_superpca, "87")
getPathSVD(colon_superpca, "KEGG_ERBB_SIGNALING_PATHWAY")
getPathSVD(colon_superpca, "pathway87")

getPathSVD(colon_aespca, "87")
getPathSVD(colon_aespca, "KEGG_ERBB_SIGNALING_PATHWAY")
getPathSVD(colon_aespca, "pathway87")
