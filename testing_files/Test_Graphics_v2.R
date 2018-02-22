# Test Plot Building
library(tidyverse)
library(reshape2)

###  Data  ###
# Grab the data from the output of the Test_supervisedPCA_wrapper.R file,
#   section "V 4 Tests". This is the 177 tumour observations


######  Plots  ################################################################

###  Categorical Response (Censoring Indicator)  ###
# Melt the data
# We have to flip the order of the pathways so the graph looks nice
classif_melt_df <- classifTest_df %>%
  head(25) %>%
  select(-terms, - setsize) %>%
  melt(id.vars = "pathways") %>%
  mutate(score = -log(value)) %>%
  mutate(pathways = factor(pathways,
                           levels = rev(unique(pathways)),
                           ordered = TRUE))

# Plot Melted Data
ggplot(classif_melt_df) +
  aes(x = pathways, y = score, fill = variable) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "FDR Adjustment",
                      breaks = c("rawp", "BH", "SidakSS"),
                      labels = c("None", "BH", "Sidak SS")) +
  ggtitle("Top 25 Tumour Pathways by Censoring Indicator") +
  xlab("Pathways") +
  ylab("Log p-Value") +
  geom_hline(yintercept = -log(0.05), size = 2) +
  geom_hline(yintercept = -log(0.1)) +
  coord_flip()



###  Regression Response (Event Time)  ###
# Melt the data
reg_melt_df <- regTest_df %>%
  head(25) %>%
  select(-terms, - setsize) %>%
  melt(id.vars = "pathways") %>%
  mutate(score = -log(value)) %>%
  mutate(pathways = factor(pathways,
                           levels = rev(unique(pathways)),
                           ordered = TRUE))

# Plot Melted Data
ggplot(reg_melt_df) +
  aes(x = pathways, y = score, fill = variable) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "FDR Adjustment",
                      breaks = c("rawp", "BH", "SidakSS"),
                      labels = c("None", "BH", "Sidak SS")) +
  ggtitle("Top 25 Tumour Pathways by Event Time") +
  xlab("Pathways") +
  ylab("Log p-Value") +
  geom_hline(yintercept = -log(0.05), size = 2) +
  geom_hline(yintercept = -log(0.1)) +
  coord_flip()



###  Full Survival Response  ###
# Melt the data
surv_melt_df <- survTest_df %>%
  head(25) %>%
  select(-terms, - setsize) %>%
  melt(id.vars = "pathways") %>%
  mutate(score = -log(value)) %>%
  mutate(pathways = factor(pathways,
                           levels = rev(unique(pathways)),
                           ordered = TRUE))

# Plot Melted Data
ggplot(surv_melt_df) +
  aes(x = pathways, y = score, fill = variable) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(name = "FDR Adjustment",
                      breaks = c("rawp", "BH", "SidakSS"),
                      labels = c("None", "BH", "Sidak SS")) +
  ggtitle("Top 25 Tumour Pathways by Survival Response") +
  xlab("Pathways") +
  ylab("Log p-Value") +
  geom_hline(yintercept = -log(0.05), size = 2) +
  geom_hline(yintercept = -log(0.1)) +
  coord_flip()
