# Reading and cleaning Alex Pico's EntrezID .gmt file.
# Gabriel Odom
# 20180601

# Alex sent us a .gmt file. I'm going to import it, clean the TERMS field, and
#   save it to the development branch of the package.

library(pathwayPCA)
path1 <- "../wikipathways-20180510-gmt-Homo_sapiens.gmt"
wikipwsHS_Entrez_pathwayCollection <- read_gmt(path1)
wikipwsTERMS_ls <- strsplit(wikipwsHS_Entrez_pathwayCollection$TERMS, split = "%")

wikipwsHS_Entrez_pathwayCollection$TERMS <- sapply(wikipwsTERMS_ls, `[[`, 3)
wikipwsHS_Entrez_pathwayCollection$description <- sapply(wikipwsTERMS_ls, `[[`, 1)
names(wikipwsHS_Entrez_pathwayCollection$pathways) <- wikipwsHS_Entrez_pathwayCollection$TERMS


# devtools::use_data(wikipwsHS_Entrez_pathwayCollection)

# I added the data and checked the package. I got this warning:
#   "Warning: found non-ASCII strings"
#   "'H19 action Rb-E2F1 signaling and CDK-N2-catenin activity' in object
#      'wikipwsHS_Entrez_pathwayCollection'"
#   "'miR-148a/miR-31/FIH1/HIF1N1-Notch signaling in glioblastoma ' in object
#      'wikipwsHS_Entrez_pathwayCollection'"
#   I did some googling, and found the "iconv()" function. I need to ask Alex
#   if this will be an issue.

# This code fails
which(wikipwsHS_Entrez_pathwayCollection$description %in%
        "H19 action Rb-E2F1 signaling and CDK-N2-catenin activity")
# This doesn't
grep("H19 action", wikipwsHS_Entrez_pathwayCollection$description)
wikipwsHS_Entrez_pathwayCollection$description[147]
grep("miR-148a", wikipwsHS_Entrez_pathwayCollection$description)
wikipwsHS_Entrez_pathwayCollection$description[353]

# This conversion fails.
iconv(wikipwsHS_Entrez_pathwayCollection$description[147],
      from = "UTF-8", to = "ASCII")
# I guess I'll try to fix it manually?
# These are the originals:
#   "H19 action Rb-E2F1 signaling and CDK-β-catenin activity"
#   "miR-148a/miR-31/FIH1/HIF1α-Notch signaling in glioblastoma "
# The greek letters are the culprits.

# My fix:
wikipwsHS_Entrez_pathwayCollection$description[147]
wikipwsHS_Entrez_pathwayCollection$description[147] <-
  "H19 action Rb-E2F1 signaling and CDK-BETA-catenin activity"

wikipwsHS_Entrez_pathwayCollection$description[353]
wikipwsHS_Entrez_pathwayCollection$description[353] <-
  "miR-148a/miR-31/FIH1/HIF1ALPHA-Notch signaling in glioblastoma "


devtools::use_data(wikipwsHS_Entrez_pathwayCollection, overwrite = TRUE)
