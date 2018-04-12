test_OS <- colon_OmicsSurv

######  Create Acessor Methods  ###############################################

###  Assay  ###
setGeneric("getAssay",
           function(object, ...){
             standardGeneric("getAssay")
           })

setMethod(f = "getAssay", signature = "OmicsPathway",
          definition = function(object, ...){
            object@assayData_df
          })

# Test
getAssay(test_OS)


###  pathwaySet  ###
setGeneric("getPathwaySet",
           function(object, ...){
             standardGeneric("getPathwaySet")
           })

setMethod(f = "getPathwaySet", signature = "OmicsPathway",
          definition = function(object, ...){
            object@pathwaySet
          })

# Test
getPathwaySet(test_OS)


###  Event Time  ###
setGeneric("getEventTime",
           function(object, ...){
             standardGeneric("getEventTime")
           })

setMethod(f = "getEventTime", signature = "OmicsSurv",
          definition = function(object, ...){
            object@eventTime
          })

# Test
getEventTime(test_OS)


###  Event Indicator  ###
setGeneric("getEvent",
           function(object, ...){
             standardGeneric("getEvent")
           })

setMethod(f = "getEvent", signature = "OmicsSurv",
          definition = function(object, ...){
            object@eventObserved
          })

# Test
getEvent(test_OS)


###  Response Vector  ###
setGeneric("getResponse",
           function(object, ...){
             standardGeneric("getResponse")
           })

setMethod(f = "getResponse", signature = "OmicsPathway",
          definition = function(object, ...){
              object@response
          })

# Test
getResponse(test_OS)
# WRITE A FULL TEST.



######  Create Replacement Methods  ###########################################
test_OS <- colon_OmicsSurv

###  Assay  ###
setGeneric("getAssay<-",
           function(object, value){
             standardGeneric("getAssay<-")
           })

setMethod(f = "getAssay<-", signature = "OmicsPathway",
          definition = function(object, value){

            object@assayData_df <- value

            if(validObject(object)){
              return(object)
            }

          })

# Test
getAssay(test_OS) <- 12
getAssay(test_OS) <- data.frame(a = 12)
getAssay(test_OS) <- data.frame(a = rep(1, 250))



###  pathwaySet  ###
setGeneric("getPathwaySet<-",
           function(object, value){
             standardGeneric("getPathwaySet<-")
           })

setMethod(f = "getPathwaySet<-", signature = "OmicsPathway",
          definition = function(object, value){

            object@pathwaySet <- value

            if(validObject(object)){
              return(object)
            }

          })

# Test
getPathwaySet(test_OS) <- 12
getPathwaySet(test_OS) <- list(a = 12)
