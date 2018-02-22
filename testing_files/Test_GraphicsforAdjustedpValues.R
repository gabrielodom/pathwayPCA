# Graphics Testing


######  AES-PCA  ######
load("inst/aespcaPathwaypVals_df.rda")

###  Data Wrangling  ###
head(aespcaPathwaypVals_df)

library(magrittr)
library(tidyverse)
aespcaPathwaypVals_df %<>% arrange(BH, rawp)
aespcaPathwaypVals_df %<>% mutate(pathScoreRaw = -log(rawp))
aespcaPathwaypVals_df %<>% mutate(pathScoreBH = -log(BH))
aespcaPathwaypVals_df$rank <- 1:nrow(aespcaPathwaypVals_df)
head(aespcaPathwaypVals_df)
tail(aespcaPathwaypVals_df)

library(ggplot2)
###  p-Values  ###
ggplot(data = aespcaPathwaypVals_df,
       aes(x = rank, y = pathScoreRaw)) +
  geom_bar(stat = "identity")
ggplot(data = aespcaPathwaypVals_df,
       aes(x = rank, y = pathScoreBH)) +
  geom_bar(stat = "identity")

###  p-Values by Pathway Size  ###
ggplot(data = aespcaPathwaypVals_df,
       aes(x = setsize, y = pathScoreRaw)) +
  geom_point()
ggplot(data = aespcaPathwaypVals_df,
       aes(x = setsize, y = pathScoreBH)) +
  geom_point()
# No relationship - good.

###  Top p-Values  ###
ggplot(data = aespcaPathwaypVals_df,
       aes(x = rank, y = pathScoreRaw, colour = BH)) +
  geom_bar(stat = "identity")
ggplot(data = aespcaPathwaypVals_df[1:500,],
       aes(x = rank, y = pathScoreRaw, colour = BH)) +
  geom_bar(stat = "identity")
ggplot(data = aespcaPathwaypVals_df[1:100,],
       aes(x = rank, y = pathScoreRaw, colour = BH)) +
  geom_bar(stat = "identity")



######  Supervised PCA  ######
data(spcaPathwayPvals_df)
head(spcaPathwayPvals_df)

###  Data Wrangling  ###
spcaPathwayPvals_df %<>% arrange(BH, rawp)
spcaPathwayPvals_df %<>% mutate(pathScoreRaw = -log(rawp))
spcaPathwayPvals_df %<>% mutate(pathScoreBH = -log(BH))
spcaPathwayPvals_df$rank <- 1:nrow(spcaPathwayPvals_df)
head(spcaPathwayPvals_df)
tail(spcaPathwayPvals_df)


###  Top p-Values  ###
ggplot(data = aespcaPathwaypVals_df,
       aes(x = rank, y = pathScoreRaw, colour = BH)) +
  geom_bar(stat = "identity")
ggplot(data = aespcaPathwaypVals_df[1:500,],
       aes(x = rank, y = pathScoreRaw, colour = BH)) +
  geom_bar(stat = "identity")
ggplot(data = aespcaPathwaypVals_df[1:100,],
       aes(x = rank, y = pathScoreRaw, colour = BH)) +
  geom_bar(stat = "identity")
