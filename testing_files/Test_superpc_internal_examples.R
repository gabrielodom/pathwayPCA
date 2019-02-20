
data("colon_pathwayCollection")
data("colonSurv_df")

colon_OmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[,-(2:3)],
  pathwayCollection_ls = colon_pathwayCollection,
  response = colonSurv_df[, 1:3],
  respType = "surv"
)

asthmaGenes_char <-
  getTrimPathwayCollection(colon_OmicsSurv)[["KEGG_ASTHMA"]]$IDs

data_ls <- list(
  x = t(getAssay(colon_OmicsSurv))[asthmaGenes_char, ],
  y = getEventTime(colon_OmicsSurv),
  censoring.status = getEvent(colon_OmicsSurv),
  featurenames = asthmaGenes_char
)

superpcFit <- superpc.train(
  data = data_ls,
  type = "surv"
)

superpc.st(
  fit = superpcFit,
  data = data_ls
)
