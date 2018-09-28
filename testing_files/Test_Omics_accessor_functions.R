
library(pathwayPCA)
data("colonSurv_df")
data("colon_pathwayCollection")

colon_OmicsPath <- CreateOmicsPath(assayData_df = colonSurv_df[, -(1:2)],
                                    pathwayCollection_ls = colon_pathwayCollection)

colon_OmicsSurv <- CreateOmics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df[, 1:2],
  respType = "survival"
)

colon_OmicsReg <- CreateOmics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df$OS_time,
  respType = "regression"
)

colon_OmicsCateg <- CreateOmics(
  assayData_df = colonSurv_df[, -(1:2)],
  pathwayCollection = colon_pathwayCollection,
  response = colonSurv_df$OS_event,
  respType = "categorical"
)

######  Create Acessor Methods  ###############################################

###  Assay  ###
# setGeneric("getAssay",
#            function(object, ...){
#              standardGeneric("getAssay")
#            })
#
# setMethod(f = "getAssay", signature = "OmicsPathway",
#           definition = function(object, ...){
#             object@assayData_df
#           })

# Test
getAssay(colon_OmicsSurv)


###  pathwayCollection  ###
# setGeneric("getPathwayCollection",
#            function(object, ...){
#              standardGeneric("getPathwayCollection")
#            })
#
# setMethod(f = "getPathwayCollection", signature = "OmicsPathway",
#           definition = function(object, ...){
#             object@pathwayCollection
#           })

# Test
getPathwayCollection(colon_OmicsSurv)


###  Event Time  ###
# setGeneric("getEventTime",
#            function(object, ...){
#              standardGeneric("getEventTime")
#            })
#
# setMethod(f = "getEventTime", signature = "OmicsSurv",
#           definition = function(object, ...){
#             object@eventTime
#           })

# Test
getEventTime(colon_OmicsPath)
getEventTime(colon_OmicsSurv)
getEventTime(colon_OmicsReg)
getEventTime(colon_OmicsCateg)


###  Event Indicator  ###
# setGeneric("getEvent",
#            function(object, ...){
#              standardGeneric("getEvent")
#            })
#
# setMethod(f = "getEvent", signature = "OmicsSurv",
#           definition = function(object, ...){
#             object@eventObserved
#           })

# Test
getEvent(colon_OmicsPath)
getEvent(colon_OmicsSurv)
getEvent(colon_OmicsReg)
getEvent(colon_OmicsCateg)


###  Response Vector  ###
# setGeneric("getResponse",
#            function(object, ...){
#              standardGeneric("getResponse")
#            })
#
# setMethod(f = "getResponse", signature = "OmicsPathway",
#           definition = function(object, ...){
#               object@response
#           })

# Test
getResponse(colon_OmicsPath)
getResponse(colon_OmicsSurv)
getResponse(colon_OmicsReg)
getResponse(colon_OmicsCateg)



######  Create Replacement Methods  ###########################################
test_OP <- colon_OmicsPath
test_OS <- colon_OmicsSurv
test_OR <- colon_OmicsReg
test_OC <- colon_OmicsCateg


###  Assay  ###
# setGeneric("getAssay<-",
#            function(object, value){
#              standardGeneric("getAssay<-")
#            })
#
# setMethod(f = "getAssay<-", signature = "OmicsPathway",
#           definition = function(object, value){
#
#             object@assayData_df <- value
#
#             if(validObject(object)){
#               return(object)
#             }
#
#           })

# Test
getAssay(test_OS) <- 12
getAssay(test_OS) <- data.frame(a = 12)
getAssay(test_OS) <- data.frame(a = rep(1, 250))


###  pathwayCollection  ###
# setGeneric("getPathwayCollection<-",
#            function(object, value){
#              standardGeneric("getPathwayCollection<-")
#            })
#
# setMethod(f = "getPathwayCollection<-", signature = "OmicsPathway",
#           definition = function(object, value){
#
#             object@pathwayCollection <- value
#
#             if(validObject(object)){
#               return(object)
#             }
#
#           })

# Test
getPathwayCollection(test_OS) <- 12
getPathwayCollection(test_OS) <- list(a = 12)
getPathwayCollection(test_OS) <-
  CreatePathwayCollection(pathways = list(a = "a0"), TERMS = "b0")


###  Event Time  ###
# setGeneric("getEventTime<-",
#            function(object, value){
#              standardGeneric("getEventTime<-")
#            })
#
# setMethod(f = "getEventTime<-", signature = "OmicsSurv",
#           definition = function(object, value){
#
#             object@eventTime <- value
#
#             if(validObject(object)){
#               return(object)
#             }
#
#           })

# Test
getEventTime(test_OS) <- 12
getEventTime(test_OS) <- rep(12, 250)


###  Event Indicator  ###
# setGeneric("getEvent<-",
#            function(object, value){
#              standardGeneric("getEvent<-")
#            })
#
# setMethod(f = "getEvent<-", signature = "OmicsSurv",
#           definition = function(object, value){
#
#             object@eventObserved <- value
#
#             if(validObject(object)){
#               return(object)
#             }
#
#           })

# Test
getEvent(test_OS) <- 12
getEvent(test_OS) <- TRUE
getEvent(test_OS) <- rep(TRUE, 250)


###  Response Vector  ###
# setGeneric("getResponse<-",
#            function(object, value){
#              standardGeneric("getResponse<-")
#            })
#
# setMethod(f = "getResponse<-", signature = "OmicsPathway",
#           definition = function(object, value){
#
#             object@response <- value
#
#             if(validObject(object)){
#               return(object)
#             }
#
#           })

# Test
getResponse(test_OP) <- 12
getResponse(test_OS) <- 12
getResponse(test_OR) <- 12
getResponse(test_OC) <- 12

getResponse(test_OR) <- rep(12, 250)

getResponse(test_OC) <- as.factor(12)
getResponse(test_OC) <- as.factor(rep(12, 250))
