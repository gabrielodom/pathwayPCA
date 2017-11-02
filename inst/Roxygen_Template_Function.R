#' Function Title
#'
#' @description A single-sentence description of the function
#'
#' @param functionArg A brief description of the argument "functionArg"
#'
#' @return What does the function return?
#'
#' @details A thorough explination of how the function works
#'
#' @export
#'   Do you want this included in the NAMESPACE file for the end-user?
#'
#' @include somefile.extension anotherfile.extension
#'   For S4 generics, we need to load the class first. The class definition is
#'   in one of those files
#'
#' @importFrom somePackage someFunction
#'   This is a function we need to call in our function, but it's in another
#'   package
#'
#' @examples egFunction(functionArg = 1)
egFunction <- function(functionArg){
  print(functionArg)
}
