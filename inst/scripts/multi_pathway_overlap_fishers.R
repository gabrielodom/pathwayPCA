# Fisher's Exact Test for Multi-Omics Pathway Overlap
# Gabriel Odom
# 2019-05-24


######  Setup  ################################################################
library(tidyverse)
library(pathwayPCA)

gitHubPath_char <-
  "https://raw.githubusercontent.com/lizhongliu1996/pathwayPCAdata/master/"
dataDir_path <- system.file(
  "extdata", package = "pathwayPCA", mustWork = TRUE
)

wikipathways_PC <- read_gmt(
  paste0(dataDir_path, "/wikipathways_human_symbol.gmt"),
  description = TRUE
)


###  Protein Expression  ###
ovSurv_df <- readRDS(
  url(paste0(gitHubPath_char, "ovarian_PNNL_survival.RDS"))
)

ov_OmicsSurv <- CreateOmics(
  assayData_df = ovSurv_df[, -(2:3)], 
  pathwayCollection_ls = wikipathways_PC,
  response = ovSurv_df[, 1:3],
  respType = "survival",
  minPathSize = 5  
)
ovProt_aespcOut <- AESPCA_pVals(
  object = ov_OmicsSurv,
  numPCs = 1,
  parallel = TRUE,
  numCores = 2,
  numReps = 0,
  adjustment = "BH"
)


###  Copy Number  ###
# copyNumberClean_df <- readRDS(
#   url(paste0(gitHubPath_char, "OV_surv_x_CNV2.RDS"))
# )
# ovCNV_Surv <- CreateOmics(
#   assayData_df = copyNumberClean_df[, -(2:3)],
#   pathwayCollection_ls = wikipathways_PC,
#   response = copyNumberClean_df[, 1:3],
#   respType = "survival",
#   minPathSize = 5
# )
# ovCNV_aespcOut <- AESPCA_pVals(
#   object = ovCNV_Surv,
#   numPCs = 1,
#   parallel = TRUE,
#   numCores = 20,
#   numReps = 0,
#   adjustment = "BH"
# )
ovCNV_aespcOut <- readRDS(
  url(paste0(gitHubPath_char, "ovarian_copyNum_aespcOut.RDS"))
)



######  Set Cardinalities  ####################################################

ovCNVpaths  <- ovCNV_aespcOut$pVals_df$description
ovProtPaths <- ovProt_aespcOut$pVals_df$description

ovCNVsignifPaths  <- getPathpVals(ovCNV_aespcOut, alpha = 0.05)$description
ovProtSignifPaths <- getPathpVals(ovProt_aespcOut, alpha = 0.05)$description


###  Total Pathways  ###
length(ovCNVpaths)  # 424
length(ovProtPaths) # 324

length(
  intersect(
    ovCNVpaths,
    ovProtPaths
  )
)
# 320 (or 324 by TERMS)

###  Significant Pathways  ###
# 1. copy number variation
length(ovCNVsignifPaths)   # 128
# 2. protein
length(ovProtSignifPaths) # 50
# 3. both protein and CNV = 22
posPos <- length(
  intersect(
    ovCNVsignifPaths,
    ovProtSignifPaths
  )
)
posPos
# (or 23 by description, because we have two WnT signalling pathways)

# 4. CNV but not protein = 106 (104)
length(
  setdiff(
    ovCNVsignifPaths,
    ovProtSignifPaths
  )
)

# 5. protein but not CNV = 28 (27)
length(
  setdiff(
    ovProtSignifPaths,
    ovCNVsignifPaths
  )
)

# 6. pathways recorded for proteins, not significant for proteins, but
#   significant for CNV = 84 (82)
negPos <- length(
  intersect(
    ovCNVsignifPaths,
    setdiff(
      ovProtPaths,
      ovProtSignifPaths
    )
  )
)
negPos

# 7. pathways recorded for proteins, not significant for CNV, but significant for
#   proteins = 28 (27)
posNeg <- length(
  intersect(
    ovProtSignifPaths,
    setdiff(
      ovCNVpaths,
      ovCNVsignifPaths
    )
  )
)
posNeg

# 8. Significant in either = 156
length(
  union(
    ovCNVsignifPaths,
    ovProtSignifPaths
  )
)
# (154 with description)

# 9. Significant in neither = 190 
length(
  intersect(
    setdiff(
      ovCNVpaths,
      ovCNVsignifPaths
    ),
    setdiff(
      ovProtPaths,
      ovProtSignifPaths
    )
  )
)
# (188 with description)


###  Fisher's Exact Test  ###
# within 324 protein pathways
# We can only look for pathways recorded in the Protein data

negNeg <- 324 - sum(c(posPos, posNeg, negPos))
negNeg

fisher.test(
  x = matrix(c(posPos, posNeg, negPos, negNeg), ncol = 2),
  alternative = "greater"
)

