# Gene Ranking by Adjusted p-Values

topFeatures <- function(object, pVals_df, percentile = 0.95){
  # browser()

  clean_obj <- expressedOmes(object, message = FALSE)
  paths_ls <- clean_obj@trimPathwaySet$pathways
  genes_char <- unique(do.call(c, paths_ls))
  rm(object, clean_obj)


  ###  Calculate Pathway Scores  ###
  ranks_df <- pVals_df[, 1, drop = FALSE]
  # Currently, the first five columns are pathways, setsize, trim_size, terms,
  #   and rawp. The adjusted p-value columns do not start until column 6.
  adj_pVals_df <- pVals_df[, 6:ncol(pVals_df), drop = FALSE]
  # Add in a buffer in case one of the p-values is identically 0 (this is
  #   probable if the number of permutations is 1000 or smaller).
  minp <- min(unlist(adj_pVals_df)[unlist(adj_pVals_df) > 0])
  ranks_df$score <- rowMeans(-log(adj_pVals_df + minp / 100))

  ###  Create a Matrix of Pathway Membership  ###
  membership_mat <- sapply(1:nrow(pVals_df), function(path){

    idx <- which(names(paths_ls) == pVals_df[path, 1])
    (genes_char %in% paths_ls[[idx]]) * ranks_df[path, 2]

  })
  rownames(membership_mat) <- genes_char


  ###  Return the Top 5% Most Significant Genes  ###
  geneScores_num <- rowSums(membership_mat)
  cutoff <- stats::quantile(geneScores_num, percentile)
  sort(geneScores_num[geneScores_num > cutoff], decreasing = TRUE)

}

# Test
# # The gene set list must come from the clean Omics* object (after it has been
# #   passed through the expressedOmes() function).
# genesetsClean_ls <- expressedOmes(colon_OmicsSurv)@pathwaySet
topFeatures(object = colon_OmicsSurv, pVals_df = surv_pVals_df)

# If you stored the results from the supervised PCA and AES-PCA runs as
#   different p-value data frames, you could find the intersecting gene set of
#   the top 5% genes from each method. Here is that intersection for the colon
#   subset:
c("PIK3CA", "PIK3R1", "PIK3CB", "PIK3R2", "PIK3R4", "PIK3C3", "MAPK1", "MAP2K1",
  "NRAS", "GRB2", "MAPK3", "HRAS", "KRAS", "SHC1", "CRK", "PIK3R5", "PIK3R3",
  "PIK3CD", "PIK3CG", "RAF1", "MAP2K2", "AKT2", "SOS1")
