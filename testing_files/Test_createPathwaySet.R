create_pathwaySet <- function(pathways, TERMS, ...){

  ###  Class Checks  ###
  if(!is.list(pathways)){
    stop("The pathways object must be a list of pathway vectors.")
  }

  if(!is.atomic(TERMS)){
    stop("The TERMS object must be an atomic vector.")
  }

  ###  Validation Checks  ###
  nPaths <- length(pathways)
  if(length(TERMS) != nPaths){
    stop("The TERMS vector must have the same length as the pathways list.")
  }

  dotNames <- names(list(...))
  if(any(c("setsize", "trim_setsize") %in% dotNames)){
    warning("The names 'setsize' and 'trim_setsize' are reserved names of a pathwaySet object.
  Values stored with these names may be overwritten during pathwayPCA function execution.
  Use with extreme caution.")
  }

  ###  Create and Return pathwaySet  ###
  out_ls <- list(pathways = pathways,
                 TERMS = TERMS,
                 ...)
  class(out_ls) <- c("pathwaySet", "list")

  out_ls

}

# Test
create_pathwaySet(pathways = 5, TERMS = "five", setsize = 1)
create_pathwaySet(pathways = list(one = 1, five = 5),
                  TERMS = list(five = "five"), setsize = 1)
create_pathwaySet(pathways = list(one = 1, five = 5),
                  TERMS = "five", setsize = 1)
create_pathwaySet(pathways = list(one = 1,
                                  two = 1:2,
                                  three = 1:3,
                                  four = 1:4,
                                  five = 1:5),
                  TERMS = c("one", "two", "three", "four", "five"),
                  setsize = 1)
# I can't do anything to check that other values of the list match the number
#   of pathways in "pathways", but I probably don't want to. If I leave it
#   flexible, we can add other stuff later.
create_pathwaySet(pathways = list(one = 1,
                                  two = 1:2,
                                  three = 1:3,
                                  four = 1:4,
                                  five = 1:5),
                  TERMS = c("one", "two", "three", "four", "five"),
                  setsize = 1:5)
