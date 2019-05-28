# biomaRt Conversion between Gene Symbol and Entrez IDs
# Gabriel Odom
# 2018-07-26



######  Create Matching Table  ################################################
# Use this walkthrough as a guide:
# https://bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html


###  Vector of Genes  ###
library(pathwayPCA)
wikipHuman_pathwaySet <- read_gmt("extdata/wikipathways_human.gmt")
allEntrezIDs_char <- unique(
  unname(unlist(wikipHuman_pathwaySet$pathways))
)
allEntrezIDs_char <- gsub("\r", "", allEntrezIDs_char)


###  BioMart Database and Dataset  ###
library(biomaRt)
# Code from James:
mart <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)


###  BioMart Query  ###
martAttr  <- listAttributes(mart = mart)
martFiltr <- listFilters(mart = mart)

humanID_df <- getBM(
  attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id", "go_id"),
  filter = "entrezgene",
  values = allEntrezIDs_char,
  mart = mart
)

write.csv(humanID_df,
          file = "gene_conversion_table.csv",
          row.names = FALSE)



######  Convert from EntrezID to Gene Symbol  #################################
###  Setup  ###
library(pathwayPCA)
wikipHuman_pathwaySet <- read_gmt("extdata/wikipathways_human.gmt")
allEntrezIDs_char <- unique(
  unname(unlist(wikipHuman_pathwaySet$pathways))
)
allEntrezIDs_char <- gsub("\r", "", allEntrezIDs_char)

library(tidyverse)
humanID_tbl <- read_csv(file = "gene_conversion_table.csv")


###  Test Conversion  ###
humanID_tbl %>%
  filter(entrezgene == as.integer(allEntrezIDs_char[1])) %>% 
  select(hgnc_symbol) %>% 
  unique %>% 
  unlist %>% 
  unname()

###  Test Vector Conversion  ###
# EntrezID 640 matches "BLK" and "NA". Other than that, we're all good
humanID_tbl %>% 
  filter(hgnc_symbol == "BLK") %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()

humanID_tbl %>% 
  filter(entrezgene == 640) %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()

humanID_tbl %>% 
  filter(entrezgene == 640) %>%
  distinct() %>% 
  as.matrix()
# The issue is that there may not be matches to the Ensembl or GO IDs. We'll
#   just omit NAs

humanID_tbl %>%
  filter(entrezgene %in% as.integer(wikipHuman_pathwaySet$pathways[[1]])) %>% 
  select(hgnc_symbol) %>% 
  na.omit %>% 
  unique %>% 
  unlist %>% 
  unname()


###  Make a Function  ###
convertGene <- function(pathway_vec,
                        conversion_tbl,
                        from_char = entrezgene,
                        to_char = hgnc_symbol){
  # browser()
  
  from_char <- enquo(from_char)
  to_char   <- enquo(to_char)
  
  conversion_tbl %>%
    filter(!!from_char %in% pathway_vec) %>% 
    select(!!to_char) %>% 
    na.omit %>% 
    unique %>% 
    unlist %>% 
    unname()
  
}

# Test
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[1]]),
  conversion_tbl = humanID_tbl
)



######  Apply the Conversion Function  ########################################

###  Setup  ###
library(pathwayPCA)
library(tidyverse)
wikipHuman_pathwaySet <- read_gmt("extdata/wikipathways_human.gmt")
humanID_tbl <- read_csv(file = "gene_conversion_table.csv")


###  Function  ###
convertGene <- function(pathway_vec,
                        conversion_tbl,
                        from_char = entrezgene,
                        to_char = hgnc_symbol){
  # browser()
  
  from_char <- enquo(from_char)
  to_char   <- enquo(to_char)
  
  conversion_tbl %>%
    filter(!!from_char %in% pathway_vec) %>% 
    select(!!to_char) %>% 
    na.omit %>% 
    unique %>% 
    unlist %>% 
    unname()
  
}


###  Apply the Function  ###
lapply(wikipHuman_pathwaySet$pathways[1:5], convertGene, humanID_tbl)

cleanWikiP_pathwaySet <- wikipHuman_pathwaySet
a <- Sys.time()
cleanWikiP_pathwaySet$pathways <- 
  lapply(wikipHuman_pathwaySet$pathways, convertGene, humanID_tbl)
Sys.time() - a # 34 sec


###  Clean Wikipathways TERMS Field  ###
terms_ls <- strsplit(wikipHuman_pathwaySet$TERMS, "%")
cleanWikiP_pathwaySet$TERMS <- sapply(terms_ls, `[`, 3)
cleanWikiP_pathwaySet$description <- sapply(terms_ls, `[`, 1)


###  Save  ###
saveRDS(cleanWikiP_pathwaySet, file = "extdata/wikipathways_human_clean.RDS")



######  Write a write_gmt Function  ###########################################

###  Setup  ###
library(pathwayPCA)
library(tidyverse)
wikipHuman_pathwaySet <- read_gmt("extdata/wikipathways_human.gmt")

test_pathwaySet <- wikipHuman_pathwaySet
test_pathwaySet$TERMS <- sapply(terms_ls, `[`, 3)
test_pathwaySet$description <- sapply(terms_ls, `[`, 1)

###  Create a File Vector  ###
# First create the list in .gmt form:
a1 <- Sys.time()
test_ls <- lapply(1:457, function(i){
  
  desc_char <- wikipHuman_pathwaySet$description[i]
  desc_char <- ifelse(is.null(desc_char), "", desc_char)
  
  c(wikipHuman_pathwaySet$TERMS[i],
    desc_char,
    wikipHuman_pathwaySet$pathways[[i]])
})
Sys.time() - a1 # 0.04464793 sec for 100; 0.04784989, 0.04785013 sec for 457

# Now collapse the list by tabs
a2 <- Sys.time()
test_char <- sapply(test_ls, paste, collapse = "\t")
Sys.time() - a2 # 0.032233 sec; 0.007007837, 0.03924084 sec for 457

# Test writing functions
a3 <- Sys.time()
writeLines(test_char, con = "extdata/test_writeLines.txt")
Sys.time() - a3 # 0.00100112 sec; 0 sec for 457

tWL <- c(0, 0.00100112, 0.03223419, 0.032233, 0.04524803)
mean(tWL)

a4 <- Sys.time()
cat(test_char, file = "extdata/test_cat.txt", sep = "\n")
Sys.time() - a4 # 0.03223395 sec; 0.03123283 sec for 457

tCat <- c(0.03123283, 0.001000881, 0.04884982, 0.06446719, 0.03123307)
mean(tCat)

# writeChar is comparable in speed to the other two methods, but doesn't write
#   the file in the proper ragged format
# a5 <- Sys.time()
# writeChar(test_char, con = "extdata/test_writeChar.txt")
# Sys.time() - a5 # x.xx sec for 457


###  Make it a Function  ###
write_gmt <- function(pathwaySet, path){
  
  ###  Setup  ###
  pathways_ls <- pathwaySet$pathways
  TERMS_char  <- pathwaySet$TERMS
  desc_char   <- pathwaySet$description
  nPaths      <- length(pathways_ls)
  
  if(is.null(desc_char)){
    desc_char <- rep(NA, nPaths)
  }
  
  ###  Write a list in .gmt form  ###
  # See the Broad Institute wiki page on .gmt file formats:
  #   https://tinyurl.com/yaf3exlk
  out_ls <- lapply(1:nPaths, function(i){
    c(TERMS_char[i], desc_char[i], pathways_ls[[i]])
  })
  
  ###  Collapse the list  ###
  out_char <- sapply(out_ls, paste, collapse = "\t")
  
  ###  Write the File  ###
  writeLines(out_char, con = path)
  
}

# Test
write_gmt(wikipHuman_pathwaySet, path = "extdata/test.gmt")



######  Write the Cleaned and Translated Wikipathways .gmt File  ##############
library(pathwayPCA)
library(tidyverse)
clean_pathwaySet <- readRDS(file = "extdata/wikipathways_human_clean.RDS")

write_gmt <- function(pathwaySet, path){
  
  ###  Setup  ###
  pathways_ls <- pathwaySet$pathways
  TERMS_char  <- pathwaySet$TERMS
  desc_char   <- pathwaySet$description
  nPaths      <- length(pathways_ls)
  
  if(is.null(desc_char)){
    desc_char <- rep(NA, nPaths)
  }
  
  ###  Write a list in .gmt form  ###
  # See the Broad Institute wiki page on .gmt file formats:
  #   https://tinyurl.com/yaf3exlk
  out_ls <- lapply(1:nPaths, function(i){
    c(TERMS_char[i], desc_char[i], pathways_ls[[i]])
  })
  
  ###  Collapse the list  ###
  out_char <- sapply(out_ls, paste, collapse = "\t")
  
  ###  Write the File  ###
  writeLines(out_char, con = path)
  
}


# Test
write_gmt(clean_pathwaySet, path = "extdata/wikipathways_human_clean.gmt")
clean2_pathwaySet <- read_gmt("extdata/wikipathways_human_clean.gmt",
                              description = TRUE)
all.equal(clean_pathwaySet, clean2_pathwaySet) # something is wrong...

all.equal(clean_pathwaySet$pathways[[1]],
          clean2_pathwaySet$pathways[[1]])
clean_pathwaySet$pathways[[1]] %in% clean2_pathwaySet$pathways[[1]]
clean_pathwaySet$pathways[[1]][96]
clean2_pathwaySet$pathways[[1]][96]
# It's the "\r" expression. It's coming from the read_gmt function, so I fixed
#   that function in the devel branch.

# Try again:
clean2_pathwaySet <- read_gmt("extdata/wikipathways_human_clean.gmt",
                              description = TRUE)
all.equal(clean_pathwaySet, clean2_pathwaySet) # less stuff is wrong, but still
all.equal(clean_pathwaySet$pathways[[167]],
          clean2_pathwaySet$pathways[[167]])
clean_pathwaySet$pathways[[167]]  # character(0)
clean2_pathwaySet$pathways[[167]] # NA and something that isn't a gene...

# What's the problem?
original_pathwaySet <- read_gmt("extdata/wikipathways_human.gmt")
original_pathwaySet$pathways[[167]] # matches gene CYP2E1 from uniprot.org
clean_pathwaySet$pathways[[167]]  # wrong
clean2_pathwaySet$pathways[[167]] # SUPER WRONG

original_pathwaySet$TERMS[[167]]
clean_pathwaySet$TERMS[[167]]
clean_pathwaySet$description[[167]]
clean2_pathwaySet$TERMS[[167]]
clean2_pathwaySet$description[[167]]

# So I need to check my lookup table for this gene.
humanID_tbl <- read_csv("gene_conversion_table.csv")
humanID_tbl %>% 
  filter(entrezgene == 1571) %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()
# So my gene table is screwy, but the write_gmt function should never write 
#   character(0) to the pathways field. Let's try it again.



######  write_gmt() Version 2  ################################################
library(pathwayPCA)
library(tidyverse)
clean_pathwaySet <- readRDS(file = "extdata/wikipathways_human_clean.RDS")

write_gmt <- function(pathwaySet, path){
  
  ###  Setup  ###
  pathways_ls <- pathwaySet$pathways
  TERMS_char  <- pathwaySet$TERMS
  desc_char   <- pathwaySet$description
  nPaths      <- length(pathways_ls)
  
  if(is.null(desc_char)){
    desc_char <- rep(NA, nPaths)
  }
  
  
  ###  Error Checking  ###
  stopifnot(nPaths == length(TERMS_char), nPaths == length(desc_char))
  pathways_ls[which(lengths(pathways_ls) == 0)] <- NA_character_
  TERMS_char[length(TERMS_char) == 0] <- NA_character_
  desc_char[length(desc_char) == 0] <- NA_character_
  
  ###  Write a list in .gmt form  ###
  # See the Broad Institute wiki page on .gmt file formats:
  #   https://tinyurl.com/yaf3exlk
  out_ls <- lapply(1:nPaths, function(i){
    c(TERMS_char[i], desc_char[i], pathways_ls[[i]])
  })
  
  ###  Collapse the list  ###
  out_char <- sapply(out_ls, paste, collapse = "\t")
  
  ###  Write the File  ###
  writeLines(out_char, con = path)
  
}

# Test
write_gmt(clean_pathwaySet, path = "extdata/wikipathways_human_clean2.gmt")
clean3_pathwaySet <- read_gmt("extdata/wikipathways_human_clean2.gmt",
                              description = TRUE)
all.equal(clean_pathwaySet, clean3_pathwaySet)
# This difference is due to the length of character(0) vs NA_character. I'm ok
#   with that.



######  Troubleshoot Matching Table  ##########################################

###  Setup  ###
library(pathwayPCA)
library(tidyverse)
wikipHuman_pathwaySet <- read_gmt("extdata/wikipathways_human.gmt")
humanID_tbl <- read_csv(file = "gene_conversion_table.csv")


###  The Issue  ###
# We have problems with the following pathways: 167, 180, 184, 244, 314, and 
#   450. They have the following entries in the original pathway set:
wikipHuman_pathwaySet$pathways[[167]] # Matching gene = CYP2E1
wikipHuman_pathwaySet$pathways[[180]] # Matching gene = CYP2E1
wikipHuman_pathwaySet$pathways[[184]] # Matching gene = CYP3A4
wikipHuman_pathwaySet$pathways[[244]] # Matching gene = MTHFR
wikipHuman_pathwaySet$pathways[[314]] # Matching gene = INS
wikipHuman_pathwaySet$pathways[[450]] # Matching gene = CYP2E1

# Can we find the correct matches?
humanID_tbl %>% 
  filter(entrezgene == 1571) %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()

humanID_tbl %>% 
  filter(entrezgene == 1576) %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()

humanID_tbl %>% 
  filter(entrezgene == 4524) %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()

humanID_tbl %>% 
  filter(entrezgene == 3630) %>% 
  select(entrezgene, hgnc_symbol) %>% 
  distinct()
# Yes.

# So what in the hell is the problem???
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[167]]),
  conversion_tbl = humanID_tbl
)
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[180]]),
  conversion_tbl = humanID_tbl
)
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[184]]),
  conversion_tbl = humanID_tbl
)
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[244]]),
  conversion_tbl = humanID_tbl
)
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[314]]),
  conversion_tbl = humanID_tbl
)
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[450]]),
  conversion_tbl = humanID_tbl
)
# My conversion function is even working properly?
convertGene(
  pathway_vec = as.integer(wikipHuman_pathwaySet$pathways[[1]]),
  conversion_tbl = humanID_tbl
)

# Rebuild "clean" list
cleanWikiP_pathwaySet <- wikipHuman_pathwaySet
cleanWikiP_pathwaySet$pathways <- 
  lapply(wikipHuman_pathwaySet$pathways, convertGene, humanID_tbl)
cleanWikiP_pathwaySet$pathways[[167]]
cleanWikiP_pathwaySet$pathways[[180]]
cleanWikiP_pathwaySet$pathways[[184]]
cleanWikiP_pathwaySet$pathways[[244]]
cleanWikiP_pathwaySet$pathways[[314]]
cleanWikiP_pathwaySet$pathways[[450]]
# these all work. What the HELL IS GOING ON???

# Clean Wikipathways TERMS Field
terms_ls <- strsplit(wikipHuman_pathwaySet$TERMS, "%")
cleanWikiP_pathwaySet$TERMS <- sapply(terms_ls, `[`, 3)
cleanWikiP_pathwaySet$description <- sapply(terms_ls, `[`, 1)

cleanWikiP2_pathwaySet <- readRDS(file = "extdata/wikipathways_human_clean.RDS")
all.equal(cleanWikiP_pathwaySet, cleanWikiP2_pathwaySet)
all.equal(cleanWikiP_pathwaySet$pathways[[1]],
          cleanWikiP2_pathwaySet$pathways[[1]])
all.equal(cleanWikiP_pathwaySet$pathways[[1]][1:89],
          cleanWikiP2_pathwaySet$pathways[[1]][1:89])
cleanWikiP_pathwaySet$pathways[[1]][90:97]
cleanWikiP2_pathwaySet$pathways[[1]][90:96]
# Where the hell did the "BCL10" gene come from?

# Also
cleanWikiP2_pathwaySet$pathways[[167]]
cleanWikiP2_pathwaySet$pathways[[180]]
cleanWikiP2_pathwaySet$pathways[[184]]
cleanWikiP2_pathwaySet$pathways[[244]]
cleanWikiP2_pathwaySet$pathways[[314]]
cleanWikiP2_pathwaySet$pathways[[450]]

# Overall, something screwy happened to these fields. I'm going to rerun the
#   write_gmt function with a browser

write_gmt <- function(pathwaySet, path){
  
  ###  Setup  ###
  pathways_ls <- pathwaySet$pathways
  TERMS_char  <- pathwaySet$TERMS
  desc_char   <- pathwaySet$description
  nPaths      <- length(pathways_ls)
  
  if(is.null(desc_char)){
    desc_char <- rep("", nPaths)
  }
  
  
  ###  Error Checking  ###
  stopifnot(nPaths == length(TERMS_char), nPaths == length(desc_char))
  pathways_ls[which(lengths(pathways_ls) == 0)] <- NA_character_
  
  ###  Write a list in .gmt form  ###
  # See the Broad Institute wiki page on .gmt file formats:
  #   https://tinyurl.com/yaf3exlk
  out_ls <- lapply(1:nPaths, function(i){
    c(TERMS_char[i], desc_char[i], pathways_ls[[i]])
  })
  
  ###  Collapse the list  ###
  out_char <- sapply(out_ls, paste, collapse = "\t")
  
  ###  Write the File  ###
  writeLines(out_char, con = path)
  
}

# Test
write_gmt(cleanWikiP_pathwaySet, path = "extdata/test.gmt")
# There's nothing wrong with this file, so I'm just going to overwrite the first
#   pathway set .gmt file I wrote.
write_gmt(cleanWikiP_pathwaySet, path = "extdata/wikipathways_human_clean.gmt")
cleanWikiP3_pathwaySet <- read_gmt("extdata/wikipathways_human_clean.gmt",
                                   description = TRUE)
all.equal(cleanWikiP_pathwaySet, cleanWikiP3_pathwaySet)
# YES!!!!!!!!!!!!!!!!!

# Ok, rename it as wikipathways_human_symbol.gmt
write_gmt(cleanWikiP_pathwaySet, path = "extdata/wikipathways_human_symbol.gmt")


# NOTE: the write_gmt function has been added to the devel branch of pathwayPCA
#   version 0.0.0.9000.