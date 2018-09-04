# Test `[[.pathwayCollection`

gmt_path <- system.file("extdata", "c2.cp.v6.0.symbols.gmt",
                        package = "pathwayPCA", mustWork = TRUE)
cp_pathwayCollection <- read_gmt(gmt_path, description = FALSE)
cp_pathwayCollection[["REACTOME_MRNA_CAPPING"]]
cp_pathwayCollection[[]]
cp_pathwayCollection[[2]]



gmt_path <- system.file("extdata", "wikipathways_human_symbol.gmt",
                        package = "pathwayPCA", mustWork = TRUE)
wikipathways_PC <- read_gmt(gmt_path, description = TRUE)
data_path <- system.file("extdata", "ovarian_PNNL_survival.RDS",
                         package = "pathwayPCA", mustWork = TRUE)
ovSurv_df <- readRDS(data_path)
ov_OmicsSurv <- create_Omics(
  assayData_df = ovSurv_df[, -(1:3)],
  pathwayCollection_ls = wikipathways_PC,
  response = ovSurv_df[, 2:3],
  respType = "survival",
  minPathSize = 5
)

wikipathways_PC[["WP75"]]
getPathwayCollection(ov_OmicsSurv)[["WP75"]]
getTrimPathwayCollection(ov_OmicsSurv)[["WP75"]]

wikipathways_PC[["WP3891"]]
getPathwayCollection(ov_OmicsSurv)[["WP3891"]]
getTrimPathwayCollection(ov_OmicsSurv)[["WP3891"]]
