# Test Parametric AESPCA p-Values with PermTestCateg
# Gabriel Odom
# 20181129

library(tidyverse)
library(pathwayPCA)

# Example data from:
# https://www.theanalysisfactor.com/generalized-linear-models-glm-r-part4/
anxiety_df <- data.frame(
  numeracy = c(
    6.6, 7.1, 7.3, 7.5, 7.9, 7.9, 8, 8.2, 8.3, 8.3, 8.4, 8.4, 8.6, 8.7, 8.8, 8.8,
    9.1, 9.1, 9.1, 9.3, 9.5, 9.8, 10.1, 10.5, 10.6, 10.6, 10.6, 10.7, 10.8, 11,
    11.1, 11.2, 11.3, 12, 12.3, 12.4, 12.8, 12.8, 12.9, 13.4, 13.5, 13.6, 13.8,
    14.2, 14.3, 14.5, 14.6, 15, 15.1, 15.7
  ),
  anxiety = c(
    13.8, 14.6, 17.4, 14.9, 13.4, 13.5, 13.8, 16.6, 13.5, 15.7, 13.6, 14, 16.1,
    10.5, 16.9, 17.4, 13.9, 15.8, 16.4, 14.7, 15, 13.3, 10.9, 12.4, 12.9, 16.6,
    16.9, 15.4, 13.1, 17.3, 13.1, 14, 17.7, 10.6, 14.7, 10.1, 11.6, 14.2, 12.1,
    13.9, 11.4, 15.1, 13, 11.3, 11.4, 10.4, 14.4, 11, 14, 13.4
  ),
  success = c(
    0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L,
    1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
  )
)

succ_mod <- glm(
  success ~ numeracy * anxiety,
  family = binomial,
  data = anxiety_df
)
succ_mod
# printing this object does not show me a p-value
str(succ_mod)
summary(succ_mod)

anova(succ_mod, test = "Chisq")
logLik(succ_mod)[1]
AIC(succ_mod)

# After reading this reminder handout on GLMs:
# http://www.stat.columbia.edu/~martin/W2024/R11.pdf,
# I need to compare this model to the null model with ANODA
null_mod <- glm(
  success ~ 1,
  family = binomial,
  data = anxiety_df
)
succ_anova <- anova(null_mod, succ_mod, test = "Chisq")
succ_anova
succ_anova[2, "Pr(>Chi)"]


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

ov_OmicsCateg <- CreateOmics(
  assayData_df = ovSurv_df[, -(1:3)],
  pathwayCollection_ls = wikipathways_PC,
  response = ovSurv_df[, 3],
  respType = "categ",
  minPathSize = 5,
  centerScale = c(TRUE, TRUE)
)


######  Test  #################################################################
# Parametric Score-based p-Value
a1 <- Sys.time()
ovarian_aespcOut <- AESPCA_pVals(
  object = ov_OmicsCateg,
  numPCs = 2,
  parallel = TRUE,
  numReps = 0,
  adjustment = "BH"
)
Sys.time() - a1 # 45.33376 sec

# Permutation-based p-Value
a2 <- Sys.time()
ovarian_aespcOut2 <- AESPCA_pVals(
  object = ov_OmicsCateg,
  numPCs = 2,
  parallel = TRUE,
  numReps = 1000,
  adjustment = "BH"
)
Sys.time() - a2 # 48.63173 sec for 100 resp (1.1x longer); 1.582723 min for
# 1000 reps (2x longer); 8.851473 min for 10k reps (12x longer)



######  Compare p-Values  #####################################################
paramP_df <-
  ovarian_aespcOut$pVals_df %>%
  select(terms, rawp, FDR_BH) %>%
  mutate(type = "score")
nonParamP_df <-
  ovarian_aespcOut2$pVals_df %>%
  select(terms, rawp, FDR_BH) %>%
  mutate(type = "permutation")

compareP_df <- inner_join(
  nonParamP_df, paramP_df, by = "terms"
)
# compareP_df <- bind_rows(nonParamP_df, paramP_df)

# Raw p-values
ggplot(data = compareP_df) +
  aes(x = -log10(rawp.x), y = -log10(rawp.y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
# scale_x_log10() +
# scale_y_log10()
cor(compareP_df$rawp.x, compareP_df$rawp.y)

# FDR
ggplot(data = compareP_df) +
  aes(x = FDR_BH.x, y = FDR_BH.y) +
  geom_point()
cor(compareP_df$FDR_BH.x, compareP_df$FDR_BH.y)

# Comment: for smaller number of permutations, the fit for FDRs is quite poor.
#   When I increased the number of reps from 100 to 1000 to 10000, the linear
#   fit actually decreased (but the correlation increased from 67% to 85%)
