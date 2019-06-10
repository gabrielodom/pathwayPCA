# Create airway Subset
# Gabriel Odom
# 2019-06-10

# The second vignette build errors on the Bioc servers.
#   Quitting from lines 190-194 (Supplement2-Importing_Data.Rmd) 
#   Error: processing vignette 'Supplement2-Importing_Data.Rmd' failed with
#   diagnostics:
#   wrong arguments for subsetting an environment
# I have added a test to make sure that this SE2Tidy() function does not break,
#   but I need to add this data subset for the vignette.

library(SummarizedExperiment)
data(airway, package = "airway")
airway_df <- SE2Tidy(airway)
airway_df <- airway_df[, 1:15]
saveRDS(airway_df, file = "vignettes/se2tidy_Exdata.RDS")
