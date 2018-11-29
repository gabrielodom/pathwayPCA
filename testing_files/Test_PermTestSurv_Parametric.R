# Test Parametric AESPCA p-Values with PermTestSurv
# Gabriel Odom
# 20181129

library(survival)
library(tidyverse)
data("lung")
str(lung)
lung_df <- lung %>%
  mutate(status = status - 1)

lung_mod <- coxph(
  Surv(time = lung_df$time, event = lung_df$status) ~ .,
  data = lung_df
)
lung_mod
# printing this object shows me a p-value, but I don't know how to access it
str(lung_mod)

summary(lung_mod)
lung_mod$loglik
lung_mod$score
lung_mod$wald.test

lungMod_summ <- summary(lung_mod)
str(lungMod_summ)
lungMod_summ$logtest
lungMod_summ$sctest
lungMod_summ$waldtest

lungMod_summ$sctest["pvalue"]

# summary(survfit(lung_mod))

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

ov_OmicsSurv <- CreateOmics(
  assayData_df = ovSurv_df[, -(1:3)],
  pathwayCollection_ls = wikipathways_PC,
  response = ovSurv_df[, 2:3],
  respType = "survival",
  minPathSize = 5,
  centerScale = c(TRUE, TRUE)
)


######  Test  #################################################################
# Parametric Score-based p-Value
a1 <- Sys.time()
ovarian_aespcOut <- AESPCA_pVals(
  object = ov_OmicsSurv,
  numPCs = 3,
  parallel = TRUE,
  numReps = 0,
  adjustment = "BH"
)
Sys.time() - a1 # 24.93709 sec
# Something's off. There are 2x as many rows of the p-values data frame (we
#   should have 457 - 133 = 324 rows). This is a problem within TabulatepValues
# EDIT: the problem wasn't in TabulatepValues, but that my edits to the
#   PermTestSurv function had the individual pathway p-value as a named atomic
#   vector. Then, when the apply call added the names to the p-value results,
#   the names were appended. Observe:
sapply(c(a = 1, b = 2), function(i) c(test = rnorm(1)))
# The two names are appended
# FIXED

# Permutation-based p-Value
a2 <- Sys.time()
ovarian_aespcOut2 <- AESPCA_pVals(
  object = ov_OmicsSurv,
  numPCs = 3,
  parallel = TRUE,
  numReps = 100,
  adjustment = "BH"
)
Sys.time() - a2 # 39.32076 sec for 100 resp (1.6x longer); 2.896029 min for
# 1000 reps (7x longer); 25.55631 min for 10k reps (61x longer)



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
  aes(x = rawp.x, y = rawp.y) +
  geom_point()
  # scale_x_log10() +
  # scale_y_log10()

# FDR
ggplot(data = compareP_df) +
  aes(x = FDR_BH.x, y = FDR_BH.y) +
  geom_point()

# Comment: for smaller number of permutations, the fit for p-values and FDRs
#   both less than 0.3 is quite poor. When I increased the number of reps from
#   1000 to 10000, the fit for (x, y) < (0.3, 0.3) improved substantially, but
#   still showed non-linearity due to the discrete nature of the non-parametric
#   p-values.
