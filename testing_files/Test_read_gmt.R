
# my.read.lines2 <- function(con){
#
#   s = file.info(con)$size
#   text_vec = readChar(con, s, useBytes = TRUE)
#   strsplit(text_vec, "\n", fixed = TRUE, useBytes = TRUE)[[1]]
#
# }



######  Speed Comparison  #####################################################
# filename <- "inst/extdata/c2.cp.v5.1.symbols.gmt"
# system.time(a_char <- scan(filename, what = character(), sep = "\t"))
# # user  system elapsed
# # 0.07    0.00    0.07
# system.time(b_char <- readLines(filename))
# # user  system elapsed
# # 0.03    0.02    0.05
# system.time(c_char <- my.read.lines2(filename))
# # user  system elapsed
# # 0.02    0.00    0.02



######  Extract First, Second, and All Other Elements  ########################
# geneset_ls <- strsplit(c_char, split = "\t")
#
# geneset_names <- sapply(geneset_ls, `[[`, 1)
# geneset_descr <- sapply(geneset_ls, `[[`, 2)
# genes_ls <- lapply(geneset_ls, function(x){
#
#   x_len <- length(x)
#   x[3:x_len]
#
# })
#
# list(genesets = genes_ls,
#      geneset.names = geneset_names,
#      geneset.descriptions = geneset_descr)



######  New Read GMT Function  ################################################
read_gmt <- function(file, delim = "\t"){

  # Read the file as a single character vector, split it by line, then split
  #   each line by
  text_char <- readChar(file, file.info(file)$size, useBytes = TRUE)
  text_vec <- strsplit(text_char, "\n", fixed = TRUE, useBytes = TRUE)[[1]]
  geneset_ls <- strsplit(text_vec, split = delim)

  # Extract the pathway names and descriptions
  geneset_names <- sapply(geneset_ls, `[[`, 1)
  geneset_descr <- sapply(geneset_ls, `[[`, 2)

  # Extract the genes
  genes_ls <- lapply(geneset_ls, function(x){

    x_len <- length(x)
    x[3:x_len]

  })

  # Create the pathwayCollection output
  out <- list(pathways = genes_ls,
              TERMS = geneset_names,
              GSEA_link = geneset_descr)
  class(out) <- "pathwayCollection"
  out

}

geneset_ls <- read_gmt("inst/extdata/c2.cp.v6.0.symbols.gmt")
system.time(geneset_ls <- read_gmt("inst/extdata/c2.cp.v6.0.symbols.gmt"))
# user  system elapsed
# 0.05    0.00    0.06
# Who's your daddy, GSA::?


######  Original (Slow) Function  #############################################
GSA.read.gmt <- function (filename){

  a = scan(filename, what = list("", ""), sep = "\t", quote = NULL,
           fill = T, flush = T, multi.line = F)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(filename, what = "", sep = "\t", quote = NULL)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
    cat(i)
    while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] !=
                                           geneset.descriptions[i])) {
      ii = ii + 1
    }
    ox[i] = ii
    ii = ii + 1
  }
  genesets = vector("list", nn)
  for (i in 1:(nn - 1)) {
    cat(i, fill = T)
    i1 = ox[i] + 2
    i2 = ox[i + 1] - 1
    geneset.descriptions[i] = dd[ox[i] + 1]
    genesets[[i]] = dd[i1:i2]
  }
  geneset.descriptions[nn] = dd[ox[nn] + 1]
  genesets[[nn]] = dd[(ox[nn] + 2):n]
  out = list(genesets = genesets, geneset.names = geneset.names,
             geneset.descriptions = geneset.descriptions)
  class(out) = "GSA.genesets"
  return(out)

}
# <environment: namespace:GSA>

genesetGSA_ls <- GSA.read.gmt("inst/extdata/c2.cp.v6.0.symbols.gmt")
system.time(genesetGSA_ls <- GSA.read.gmt("inst/extdata/c2.cp.v6.0.symbols.gmt"))
# user  system elapsed
# 1.08    0.05    1.09
