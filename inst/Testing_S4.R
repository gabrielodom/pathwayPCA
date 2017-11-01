

######  Create S4 Classes  ####################################################

library(methods)
load("data/ovarianFiltered_df.rda")
load("data/genesets_ls.rda")

# Create a Master Class for all experiments

# so my create.*() functions don't throw a fit when they see a tibble:
setOldClass(c("tbl_df", "tbl", "data.frame"))

create_OmicsPath <- setClass("OmicsPathway",
                             slots = c(massSpec = "data.frame",
                                       pathwaySet = "list"))

###  Now initialise creator functions for each of the SRC data types  ###

# Survival (numeric and factor response)
create_OmicsSurv <- setClass("OmicsSurv",
                             slots = c(eventTime = "numeric",
                                       eventObserved = "logical"),
                             contains = "OmicsPathway")
Y_time <- rnorm(58, mean = 78, sd = 6)
Y_event <- sample(c(FALSE, TRUE), 58, replace = TRUE, prob = c(0.2, 1 - 0.2))
testOmicsSurv <- create_OmicsSurv(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = pathways,
                                  eventTime = Y_time,
                                  eventObserved = Y_event)
# getClass(testOmicsSurv)

# Regression (continuous response)
create_OmicsReg <- setClass("OmicsReg",
                            slots = c(response = "numeric"),
                            contains = "OmicsPathway")
Y_reg <- rnorm(58)
testOmicsReg <- create_OmicsReg(massSpec = ovarianFiltered_df[, -(1:3)],
                                pathwaySet = pathways,
                                response = Y_reg)

# Classification (factor response)
create_OmicsClassif <- setClass("OmicsClassif",
                                slots = c(response = "factor"),
                                contains = "OmicsPathway")
Y_class <- factor(sample(c("A", "B", "C"), 58, replace = TRUE))
testOmicsClassif <- create_OmicsClassif(massSpec = ovarianFiltered_df[, -(1:3)],
                                        pathwaySet = pathways,
                                        response = Y_class)

# Pathway Extraction only
testOmicsPath <- create_OmicsPath(massSpec = ovarianFiltered_df[, -(1:3)],
                                  pathwaySet = pathways)


######  Create Functions for these S4 Classes  ################################

###  Geneset Extraction Function  ###
setGeneric("expressedOmes",
           function(object, ...){
             standardGeneric("expressedOmes")
             }
           )
setMethod(f = "expressedOmes", signature = "OmicsPathway",
          definition = function(object, trim = 3){
            genelist <- colnames(object@massSpec)

            # Delete the genes from the pathways if they aren't recorded in our
            #   data matrix
            newPaths <- sapply(object@pathwaySet, function(x){
              x[x %in% genelist]
            })

            # Remove any pathway that now has fewer than "trim" genes
            newPaths_trim <- sapply(seq_along(newPaths), function(i){

              shortPath <- length(newPaths[[i]]) < trim
              if(shortPath){
                NULL
              } else {
                newPaths[[i]]
              }

            })

            names(newPaths_trim) <- names(object@pathwaySet)

            nullPaths <- which(sapply(newPaths_trim, is.null))
            paths_ls <- newPaths_trim[-nullPaths]
            attr(paths_ls, "missingPaths") <- names(newPaths_trim)[nullPaths]
            paths_ls

          })

# Test:
testRedPath <- expressedOmes(testOmicsPath)
testRedPath2 <- expressedOmes(testOmicsSurv)
testRedPath3 <- expressedOmes(testOmicsReg)
testRedPath4 <- expressedOmes(testOmicsClassif)
identical(testRedPath, testRedPath2)
identical(testRedPath, testRedPath3)
identical(testRedPath, testRedPath4)
# It works
