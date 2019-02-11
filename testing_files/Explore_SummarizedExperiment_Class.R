# Explore SummarizedExperiment data class
# Gabriel Odom
# 20190208

library(SummarizedExperiment)
data(airway, package = "airway")

######  Explore the Class  ######
str(airway)
slotNames(airway)
# "rowRanges" "colData" "assays" "NAMES" "elementMetadata" "metadata"

###  rowRanges  ###
airway@rowRanges
# This is a list of GRanges (data frame-ish?) objects: one for each feature in
#   the data (64102 of them). Features are columns in this package.
str(airway@rowRanges)
# This isn't even a normal list, it's a GRangesList
slotNames(airway@rowRanges)
# "unlistData" "elementMetadata" "elementType" "metadata" "partitioning"
airway@rowRanges@unlistData
# 745593 "ranges" with 2 metadata columns each
slotNames(airway@rowRanges@unlistData)
# "seqnames" "ranges" "strand" "seqinfo" "elementMetadata" "elementType"
#   "metadata"

# I think I need to focus on supporting SummarizedExperiment, rather than its
#   "Ranged" child class.

# Here's how the rowData function works:
rowData(airway)
showMethods("rowData")
getMethod("rowData", "SummarizedExperiment")
# function (x, use.names = TRUE, ...)
#   mcols(x, use.names = use.names, ...)
# That's it? one line call to S4Vectors::mcols()?


###  colData  ###
airway@colData
# This is sample metadata
colData(airway)
showMethods("colData")
getMethod("colData", "SummarizedExperiment")
# function (x, ...)
#   x@colData


###  assays  ###
airway@assays
# R6 object of class "ShallowSimpleListAssays" with Field "data". The only slot
#   is an environment named ".xData".
ls(airway@assays@.xData)
# "data" "field" "show"
mget("data", airway@assays@.xData)$data[[1]]

assay(airway)
showMethods("assay")
getMethod("assay", c("SummarizedExperiment", "character"))
getMethod("assay", c("SummarizedExperiment", "missing"))
getMethod("assay", c("SummarizedExperiment", "numeric"))
# We don't plan to subset the data here, so we only need the "missing" method:
# function (x, i, ...)
# {
#   assays <- assays(x, ...)
#   if (0L == length(assays))
#     stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
#          "length(assays(<", class(x), ">)) is 0'")
#   assays[[1]]
# }

showMethods("assays")
getMethod("assays", "SummarizedExperiment")
# function (x, ..., withDimnames = TRUE)
# {
#   assays <- as(x@assays, "SimpleList")
#   if (withDimnames) {
#     assays <- endoapply(assays, function(assay) {
#       dimnames(assay)[1:2] <- dimnames(x)
#       assay
#     })
#   }
#   assays
# }

# Here, we need access to the "SimpleList" class (from S4Vectors) and
#   S4Vectors::endoapply().



######  Create New SummarizedExperiment  ######
# Sample info
nSamps <- 50
sampID <- paste0("subj_", seq_len(nSamps))
sampInfo <- data.frame(
  Treatment = rep(c("Case", "Ctrl"), each = 25),
  row.names = sampID
)

# Gene info
nGenes <- 1000
geneID <- sapply(seq_len(nGenes), function(i){

  alphas <- paste(sample(LETTERS, size = 4, replace = TRUE), collapse = "")
  nums <- round(runif(1, min = 0.5, max = 9.5))

  paste0(alphas, nums)

})
geneInfo <- data.frame(
  featureIDs = geneID,
  row.names = geneID
)

# Expression
expr_mat <- matrix(runif(n = nGenes * nSamps), ncol = nSamps)

# Data Container
test_SE <- SummarizedExperiment(
  assays = list(expression = expr_mat),
  rowData = geneInfo,
  colData = sampInfo
)



######  Transform SE -> Omics  ######
# We need the "SimpleList" class and the endoapply() function from S4Vectors
#   Just kidding, we don't need either of those?
SE2Tidy <- function(summExperiment, whichAssay = 1){
  # browser()

  ###  Assay  ###
  assay_mat <- summExperiment@assays[[whichAssay]]
  dimnames(assay_mat)[1:2] <- dimnames(summExperiment)
  assay_df <- pathwayPCA::TransposeAssay(
    assay_df = data.frame(assay_mat),
    omeNames = "rowNames"
  )

  ###  Response  ###
  samps_df <- as.data.frame(summExperiment@colData)

  ###  Return Joined Data Frame  ###
  out_df <- merge(samps_df, assay_df, by = "row.names")
  class(out_df) <- c("tbl_df", "tbl", "data.frame")

  out_df

}
airway_df <- SE2Tidy(airway)
test_df <- SE2Tidy(test_SE)
