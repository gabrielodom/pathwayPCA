
######  Assay Data  ###########################################################
mu_ls <- c(5.471388889, 9.143888889, 7.7, 5.484722222, 8.196666667, 8.5975,
           6.575, 4.902777778, 7.922777778, 6.484444444, 10.20222222,
           9.688055556, 7.378055556, 9.993333333, 9.839166667, 7.014444444,
           7.66)

sds_ls <- c(0.241189348, 0.858254966, 0.603390421, 0.400481754, 0.541859233,
            0.388450217, 0.170654203, 0.146728342, 0.212691744, 0.77408174,
            0.355649194, 0.356601574, 0.677600048, 0.374135193, 0.606019212,
            0.531646738, 0.228172867)

samps_ls <- purrr::map2(mu_ls, sds_ls, rnorm, n = 36)
names(samps_ls) <- c("SOAT1", "LSS", "SQLE", "EBP", "CYP51A1", "DHCR7",
                     "CYP27B1", "DHCR24", "HSD17B7", "MSMO1", "FDFT1", "SC5DL",
                     "LIPA", "CEL", "TM7SF2", "NSDHL", "SOAT2")
samps_mat <- t(as.matrix(dplyr::bind_cols(samps_ls)))
colnames(samps_mat) <- paste0("T211013", 11:46)

write.csv(samps_mat, file = "inst/extdata/ex_assay_subset.csv")



######  Response Data  ########################################################
assay_df <- readr::read_csv("inst/extdata/ex_assay_subset.csv")
assayT_df <- transpose_assay(assay_df)

beta <- rnorm(ncol(assayT_df) - 1, sd = 0.1)
y <- apply(as.matrix(assayT_df[, -1]), 1, function(row){
  row %*% beta + rnorm(1, sd = 1.0)
})
survMonths <- 5 * (max(y) - y) + 1
survMonths <- ceiling(survMonths * 4) / 4 # Round to nearest week
plot(survMonths, ylim = c(0, 36))

event_logi <- as.logical(runif(36) < 0.85)
plot(event_logi)

pInfo_df <- data.frame(Sample = paste0("T211013", 11:46),
                       eventTime = survMonths,
                       eventObserved = event_logi,
                       stringsAsFactors = FALSE)

write.csv(pInfo_df, file = "inst/extdata/ex_pInfo_subset.csv")



######  Test Joined Data  #####################################################
assay_df <- readr::read_csv("inst/extdata/ex_assay_subset.csv")
assayT_df <- transpose_assay(assay_df)
pInfo_df <- readr::read_csv("inst/extdata/ex_pInfo_subset.csv")

test_df <- full_join(pInfo_df, assayT_df, by = "Sample")
library(survival)
coxph(Surv(eventTime, eventObserved) ~ ., data = test_df[, -1])
