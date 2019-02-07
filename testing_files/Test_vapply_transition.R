# Testing vapply() v sapply()
# 20190207
# Gabriel Odom


data("colonSurv_df")
data("colon_pathwayCollection")

###  Create an OmicsSurv Object  ###
colon_Omics <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, 1:3],
  respType = "surv"
)

######  Integer / Numeric  ######
colonPCs_ls <- ExtractAESPCs(
  object = colon_Omics,
  numPCs = 2,
  parallel = TRUE,
  numCores = 2
)

PermTestSurv(
  OmicsSurv = colon_Omics,
  pathwayPCs_ls = colonPCs_ls$PCs
)

pValues1_vec <- sapply(
  pathwayPCs_ls,
  permute_SurvFit,
  resp_Surv = response,
  numReps_int = numReps
)

pValues2_vec <- vapply(
  pathwayPCs_ls,
  permute_SurvFit,
  resp_Surv = response,
  numReps_int = numReps,
  FUN.VALUE = numeric(1)
)

identical(pValues1_vec, pValues2_vec)


tScoreMax_num <- sapply(tScores_ls, function(ls) absMax(ls$tscor[1, ]) )
tScoreMax_num2 <- vapply(
  tScores_ls,
  function(ls) absMax(ls$tscor[1, ]) ,
  FUN.VALUE = numeric(1)
)
identical(tScoreMax_num, tScoreMax_num2)

tControlMax_num <- sapply(tControl_ls, function(x) absMax(x[1, ]) )
tControlMax_num2 <- vapply(
  tControl_ls,
  function(x) absMax(x[1, ]) ,
  FUN.VALUE = numeric(1)
)
identical(tControlMax_num, tControlMax_num2)


######  Character  #######
###  Pathway p-Values  ###
pVals <- PermTestSurv(
  OmicsSurv = colon_Omics,
  pathwayPCs_ls = colonPCs_ls$PCs,
  parallel = TRUE,
  numCores = 2
)

###  Create Table of p-Values  ###
trimmed_PC <- getTrimPathwayCollection(colon_Omics)
TabulatepValues(
  pVals_vec = pVals,
  genesets_ls = trimmed_PC
)
newAdjNames <- sapply(orderedNames, function(i){
  ifelse(i %in% fwerNames,
         paste0("FWER_", i),
         paste0("FDR_", i))
})
newAdjNames2 <- vapply(
  orderedNames,
  function(i){
    ifelse(
      i %in% fwerNames,
      paste0("FWER_", i),
      paste0("FDR_", i)
    )
  },
  FUN.VALUE = character(1)
)

identical(newAdjNames, newAdjNames2)


# so, as long as you specify the class of the output, vapply() is ok. I wonder
#   how it handles matrix / data frame output?
