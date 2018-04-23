######  Remove all non-alphanumeric  ##########################################

gsub("[^a-zA-Z/d]", "", "summer")
gsub("[^a-zA-Z/d]", "", "summer69")
gsub("[^a-zA-Z0-9]", "", "summer69")
gsub("[^a-zA-Z0-9]", "", "summer69.")
gsub("[^a-zA-Z0-9]", "", "summer69+")
gsub("[^a-zA-Z0-9]", "", "summer69-")
gsub("[^a-zA-Z0-9]", "", "summer69*")
gsub("[^a-zA-Z0-9]", "", "summer69/")
gsub("[^a-zA-Z0-9]", "", "summer69`")
gsub("[^a-zA-Z0-9]", "", "summer69~")
gsub("[^a-zA-Z0-9]", "", "summer69!")
gsub("[^a-zA-Z0-9]", "", "summer69@")
gsub("[^a-zA-Z0-9]", "", "summer69#")
gsub("[^a-zA-Z0-9]", "", "summer69$")
gsub("[^a-zA-Z0-9]", "", "summer69^")
gsub("[^a-zA-Z0-9]", "", "summer69&")
gsub("[^a-zA-Z0-9]", "", "summer69(")
gsub("[^a-zA-Z0-9]", "", "summer69)")
gsub("[^a-zA-Z0-9]", "", "summer69_")
gsub("[^a-zA-Z0-9]", "", "summer69=")
gsub("[^a-zA-Z0-9]", "", "summer69<")
gsub("[^a-zA-Z0-9]", "", "summer69>")
gsub("[^a-zA-Z0-9]", "", "summer69,")
gsub("[^a-zA-Z0-9]", "", "summer69?")
gsub("[^a-zA-Z0-9]", "", "summer69:")
gsub("[^a-zA-Z0-9]", "", "summer69;")
gsub("[^a-zA-Z0-9]", "", "summer69'")
gsub("[^a-zA-Z0-9]", "", "summer69"")  # Error "
gsub("[^a-zA-Z0-9]", "", "summer69{")
gsub("[^a-zA-Z0-9]", "", "summer69[")
gsub("[^a-zA-Z0-9]", "", "summer69}")
gsub("[^a-zA-Z0-9]", "", "summer69]")
gsub("[^a-zA-Z0-9]", "", "summer69|")
gsub("[^a-zA-Z0-9]", "", "summer69\")  # Error \



######  Remove Quote  #########################################################

gsub(""", "", "summer69"")
gsub('"', '', 'summer69"')   # "

test_char <- 'summer69"'
gsub('"', '', test_char)

gsub("'", "", "summer'69")



######  Remove Backslash  #####################################################

test2_char <- 'summer\69'
gsub('\', '', test2_char)     # breaks
gsub('\\', '', test2_char)    # Error
gsub('\\\', '', test2_char)   # breaks
gsub('\\\\', '', test2_char)  # doesn't remove the "\"

# I spoke with James: he said not to worry about trying to escape the "\". "


######  String Clean Function  ################################################
clean_names <- function(string){

  # WARNING: THIS FUNCTION CANNOT REMOVE "\".

  string1 <- gsub("'", "", string)  # remove any unmatched single quotes
  string2 <- gsub('"', '', string1)  # remove any unmatched double quotes
  string3 <- gsub("[^a-zA-Z0-9]", "", string2)  # remove any non-alphanumerics
  string3

}
clean_names('summer"69')
clean_names('~`!@#$%^&*()_+-={}|[]:;<>?,./summer"69')
clean_names("summer'69")

# So, I can remove a single unmatch double quote OR a single unmatched single
#   quote, but not both at the same time.


######  Detect Leading Numeral  ###############################################
test3_char <- "1Summerinthe60s"
grepl("^\\d", test3_char)
grepl("^\\d", test2_char)


######  Detect Problem Strings  ###############################################
detect_invalid_names <- function(string_ls){
  # browser()

  unmatchedSingleQ <- grepl("'", string_ls)
  unmatchedDoubleQ <- grepl('"', string_ls)
  nonAlphaNum_idx <- grepl("[^a-zA-Z0-9]", string_ls)
  leadingNum_idx <- grepl("^\\d", string_ls)

  bad_idx <- unmatchedSingleQ +
    unmatchedDoubleQ +
    nonAlphaNum_idx +
    leadingNum_idx
  bad_idx <- as.logical(bad_idx)

  string_ls[bad_idx]

}

# Test
testChar_ls <- list(a = "summer69", b = "summer69+",
                    c = 'summer"69', d = "1a",
                    e = "`yes`")
detect_invalid_names(test_char)
detect_invalid_names(testChar_ls)
detect_invalid_names(unlist(testChar_ls))



######  Include Dashes  #######################################################
# After inspecting some of my real data, I believe that proper gene names can
#   also include dashes.

testChar2_ls <- list(a = "HLA-DRB5", b = "PAK3")
grepl("[^a-zA-Z0-9]", testChar2_ls)
grepl("[^-a-zA-Z0-9-]", testChar2_ls)
