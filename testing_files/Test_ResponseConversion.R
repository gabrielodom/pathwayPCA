# Test Input Conversion

######  Regression  ###########################################################
convertRegResponse <- function(object){
  # browser()

  if(is.data.frame(object)){
    object <- as.matrix(object)
  }

  if(is.matrix(object)){

    if(ncol(object) == 1 || nrow(object) == 1){
      object <- as.vector(object)
    } else {
      stop("Multivariate response not supported at this time.")
    }

  }

  if(is.list(object)){

    if(length(object) == 1){
      object <- unlist(object, use.names = FALSE)
    } else {
      stop("Regression response must be a vector or 1-dimensional data frame.")
    }
  }

  if(is.atomic(object)){
    object <- as.numeric(object)
  } else {
    stop("Regression response must be a vector or 1-dimensional data frame.")
  }

  object

}


# Test
convertRegResponse(object = joinedExperiment_df[, 2])
convertRegResponse(object = as.matrix(joinedExperiment_df[, 2]))
convertRegResponse(object = list(time = joinedExperiment_df$eventTime))
convertRegResponse(object = list(time = joinedExperiment_df$eventTime,
                                 observed = "a"))
convertRegResponse(object = joinedExperiment_df$eventTime)



######  Classification  #######################################################
convertCatResponse <- function(object){
  # browser()

  if(is.data.frame(object)){
    object <- as.matrix(object)
  }

  if(is.matrix(object)){

    if(ncol(object) == 1 || nrow(object) == 1){
      object <- as.vector(object)
    } else {
      stop("Multivariate response not supported at this time.")
    }

  }

  if(is.list(object)){

    if(length(object) == 1){
      object <- unlist(object, use.names = FALSE)
    } else {
      stop("Regression response must be a vector or 1-dimensional data frame.")
    }
  }

  if(is.atomic(object)){
    object <- as.factor(object)
  } else {
    stop("Regression response must be a vector or 1-dimensional data frame.")
  }

  object

}


# Test
convertCatResponse(object = joinedExperiment_df[, 3])
convertCatResponse(object = as.matrix(joinedExperiment_df[, 3]))
convertCatResponse(object = joinedExperiment_df$eventObserved)
convertCatResponse(object = list(observed = joinedExperiment_df$eventObserved))
convertCatResponse(object = list(observed = joinedExperiment_df$eventObserved,
                                 time = "a"))
convertCatResponse(object = as.factor(joinedExperiment_df$eventObserved))



######  Survival  #############################################################
convertSurvResponse <- function(object){
  # browser()

  if(is.Surv(object)){
    object <- as.matrix(object)
  }

  if(is.data.frame(object)){
    object <- as.matrix(object)
  }

  if(is.matrix(object)){

    if(ncol(object) == 2){

      outTime <- object[, 1, drop = TRUE]
      outDeath <- object[, 2, drop = TRUE]

    } else {
      stop("Object must have two columns only: death time and death indicator.")
    }

  } else if(is.list(object)){

    if(length(object) == 2){

      outTime <- object[[1]]
      outDeath <- object[[2]]

    } else {
      stop("Object must have two entries only: death time and death indicator.")
    }
  }

  if(is.atomic(outTime) && is.atomic(outDeath)){

    outTime <- as.numeric(outTime)
    outDeath <- as.logical(outDeath)

  } else {
    stop("Death time and death indicator must be atomic vectors.")
  }

  list(eventTime = outTime,
       eventObserved = outDeath)

}


# Test
library(survival)
test_surv <- Surv(joinedExperiment_df$eventTime,
                  joinedExperiment_df$eventObserved)
is.Surv(test_surv)
test_surv[1, ]
test_surv[1, 2]

convertSurvResponse(test_surv)
convertSurvResponse(joinedExperiment_df[, 2:3])
convertSurvResponse(joinedExperiment_df[, 2:4])
convertSurvResponse(object = list(time = joinedExperiment_df$eventTime,
                                  observed = joinedExperiment_df$eventObserved))
convertSurvResponse(object = list(time = joinedExperiment_df$eventTime,
                                  time2 = joinedExperiment_df$eventTime,
                                  observed = joinedExperiment_df$eventObserved))



######  Combined  #############################################################

convertResponse <- function(object,
                            type = c("survival", "regression", "categorical")){

  type <- match.arg(type)

  if(type == "survival"){

    if(is.Surv(object)){
      object <- as.matrix(object)
    }

    if(is.data.frame(object)){
      object <- as.matrix(object)
    }

    if(is.matrix(object)){

      if(ncol(object) == 2){

        outTime <- object[, 1, drop = TRUE]
        outDeath <- object[, 2, drop = TRUE]

      } else {
        stop("Object must have two columns only: death time and death indicator.")
      }

    } else if(is.list(object)){

      if(length(object) == 2){

        outTime <- object[[1]]
        outDeath <- object[[2]]

      } else {
        stop("Object must have two entries only: death time and death indicator.")
      }
    }

    if(is.atomic(outTime) && is.atomic(outDeath)){

      outTime <- as.numeric(outTime)
      outDeath <- as.logical(outDeath)

    } else {
      stop("Death time and death indicator must be atomic vectors.")
    }

    list(time = outTime,
         dead = outDeath)

  } else {

    if(is.data.frame(object)){
      object <- as.matrix(object)
    }

    if(is.matrix(object)){

      if(ncol(object) == 1 || nrow(object) == 1){
        object <- as.vector(object)
      } else {
        stop("Multivariate response not supported at this time.")
      }

    }

    if(is.list(object)){

      if(length(object) == 1){
        object <- unlist(object, use.names = FALSE)
      } else {
        stop("Non-survival response must be a vector or 1-dimensional data frame.")
      }
    }

    if(is.atomic(object)){

      if(type == "regression"){
        object <- as.numeric(object)
      } else {
        object <- as.factor(object)
      }

    } else {
      stop("Non-survival response must be a vector or 1-dimensional data frame.")
    }

    object

  }

}

# Test
convertResponse(object = joinedExperiment_df[, 2], type = "r")
convertResponse(object = as.matrix(joinedExperiment_df[, 2]), type = "r")
convertResponse(object = list(time = joinedExperiment_df$eventTime), type = "r")
convertResponse(object = list(time = joinedExperiment_df$eventTime,
                                 observed = "a"),
                   type = "r")
convertResponse(object = joinedExperiment_df$eventTime, type = "r")


convertResponse(object = joinedExperiment_df[, 3], type = "c")
convertResponse(object = as.matrix(joinedExperiment_df[, 3]), type = "c")
convertResponse(object = joinedExperiment_df$eventObserved, type = "c")
convertResponse(object = list(observed = joinedExperiment_df$eventObserved),
                type = "c")
convertResponse(object = list(observed = joinedExperiment_df$eventObserved,
                                 time = "a"),
                type = "c")
convertResponse(object = as.factor(joinedExperiment_df$eventObserved),
                type = "c")


library(survival)
test_surv <- Surv(joinedExperiment_df$eventTime,
                  joinedExperiment_df$eventObserved)
convertResponse(test_surv, type = "s")
convertResponse(joinedExperiment_df[, 2:3], type = "s")
convertResponse(joinedExperiment_df[, 2:4], type = "s")
convertResponse(object = list(time = joinedExperiment_df$eventTime,
                              observed = joinedExperiment_df$eventObserved),
                type = "s")
convertResponse(object = list(time = joinedExperiment_df$eventTime,
                              time2 = joinedExperiment_df$eventTime,
                              observed = joinedExperiment_df$eventObserved),
                type = "s")
