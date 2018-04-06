######  Print method for OmicsPath objects  ###################################

# setGeneric(name = "print",
#            def = function(object, ...){
#              standardGeneric("print")
#            },
#            useAsDefault = function(object, ...){
#              base::print(object, ...)
#            })
#
# setMethod(f = "print", signature = "OmicsPathway",
#           definition = function(object, ...){
#             str(object, max.level = 2, ...)
#           })
#
# # Test
# print(test_OmicsSurv)
# # print() dispatches to show() for S4 objects


######  Show method for OmicsPath objects  ####################################
# setGeneric(name = "show",
#            def = function(object){
#              standardGeneric("show")
#            },
#            useAsDefault = function(object){
#              methods::show(object)
#            })

setMethod(f = "show", signature = "OmicsSurv",
          definition = function(object){
            str(object, max.level = 2)
          })

# Test
show(test_OmicsSurv)
# works; can't use ... in the arguments of the function, because show() doesn't
#   take ...
test_OmicsSurv
# ~~Doesn't work~~ Now it does. Don't define the generic and it works.

print(test_OmicsSurv)
# Even works for print()


### Try it with the signature  = "OmicsPathway"  ###
setMethod(f = "show", signature = "OmicsPathway",
          definition = function(object){
            str(object, max.level = 2)
          })

# Test
show(test_OmicsSurv)
test_OmicsSurv
print(test_OmicsSurv)
# It works.



######  Print for pathwaySet objects  #########################################
print <- function(x) UseMethod("print")

print.pathwaySet <- function(x){
  cat("Informal Class 'pathwaySet' [package 'pathwayPCA'] with",
      length(x), "elements: \n")
  str(x,
      max.level = 1,
      vec.len = 1,
      give.attr = FALSE,
      no.list = TRUE)
}

print(gene_set_ls)

# print.OmicsPathway <- function(x){
#   str(x, max.level = 2)
# }
# print.OmicsSurv <- print.OmicsPathway
# print.OmicsReg <- print.OmicsPathway
# print.OmicsCateg <- print.OmicsPathway
