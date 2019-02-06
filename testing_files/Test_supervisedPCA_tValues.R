# Test PathwaytValues
# 20190206
# Gabriel Odom

data("colon_pathwayCollection")
data("colonSurv_df")


######  Survival  #############################################################
colon_OmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, 1:3],
  respType = "surv"
)

asthmaGenes_char <-
  getTrimPathwayCollection(colon_OmicsSurv)[["KEGG_ASTHMA"]]$IDs
resp_mat <- matrix(
  c(getEventTime(colon_OmicsSurv), getEvent(colon_OmicsSurv)),
  ncol = 2
)

###  pathway_tScores  ###
outS1 <- PathwaytValues(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsSurv)),
  response_mat = resp_mat,
  responseType = "survival"
)
outS2 <- pathway_tScores(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsSurv)),
  response_mat = resp_mat,
  responseType = "survival"
)

all.equal(outS1, outS2)



###  pathway_tControl  ###
set.seed(1)
outC1 <- PathwaytValues(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsSurv)),
  response_mat = resp_mat,
  responseType = "survival",
  control = TRUE
)
set.seed(1)
outC2 <- pathway_tControl(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsSurv)),
  response_mat = resp_mat,
  responseType = "survival"
)

all.equal(outC1, outC2)



######  Regression  ###########################################################
colon_OmicsReg <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, 1:2],
  respType = "reg"
)

resp_mat <- matrix(
  getResponse(colon_OmicsReg),
  ncol = 1
)

###  pathway_tScores  ###
outS1 <- PathwaytValues(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsReg)),
  response_mat = resp_mat,
  responseType = "regression"
)
outS2 <- pathway_tScores(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsReg)),
  response_mat = resp_mat,
  responseType = "regression"
)

all.equal(outS1, outS2)



###  pathway_tControl  ###
set.seed(1)
outC1 <- PathwaytValues(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsReg)),
  response_mat = resp_mat,
  responseType = "regression",
  control = TRUE
)
set.seed(1)
outC2 <- pathway_tControl(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsReg)),
  response_mat = resp_mat,
  responseType = "regression"
)

all.equal(outC1, outC2)



######  Categorical  ##########################################################
colon_OmicsCateg <- CreateOmics(
  assayData_df = colonSurv_df[, -(2:3)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, c(1, 3)],
  respType = "categ"
)

resp_mat <- getResponse(colon_OmicsCateg)
dim(resp_mat) <- c(250, 1)

###  pathway_tScores  ###
outS1 <- PathwaytValues(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsCateg)),
  response_mat = resp_mat,
  responseType = "categorical"
)
outS2 <- pathway_tScores(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsCateg)),
  response_mat = resp_mat,
  responseType = "categorical"
)

all.equal(outS1, outS2)


###  pathway_tControl  ###
set.seed(1)
outC1 <- PathwaytValues(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsCateg)),
  response_mat = resp_mat,
  responseType = "categorical",
  control = TRUE
)
set.seed(1)
outC2 <- pathway_tControl(
  pathway_vec = asthmaGenes_char,
  geneArray_df = t(getAssay(colon_OmicsCateg)),
  response_mat = resp_mat,
  responseType = "categorical"
)

all.equal(outC1, outC2)
