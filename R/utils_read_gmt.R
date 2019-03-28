#' Read a \code{.gmt} file in as a \code{pathwayCollection} object
#'
#' @description Read a set list file in Gene Matrix Transposed (\code{.gmt})
#'    format, with special performance consideration for large files. Present
#'    this object as a \code{pathwayCollection} object.
#'
#' @param file A path to a file or a connection. This file must be a \code{.gmt}
#'    file, otherwise input will likely be nonsense. See the "Details" section
#'    for more information.
#' @param setType What is the type of the set: pathway set of gene, gene sites
#'    in RNA or DNA, or regions of CpGs. Defaults to \code{''pathway''}.
#' @param description Should the "description" field (the second field in the
#'    \code{.gmt} file on each line) be included in the output? Defaults to
#'    \code{FALSE}.
#' @param delim The \code{.gmt} delimiter. As proper \code{.gmt} files are tab
#'    delimited, this defaults to \code{"\\t"}.
#'
#' @return A \code{pathwayCollection} list of sets. This list has three
#'    elements:
#' \itemize{
#'   \item{'setType' : }{A named list of character vectors. Each vector
#'      contains the names of the individual genes, sites, or CpGs within that
#'      set as a vector of character strings. The name of this list entry is
#'      equal to the value specified in \code{setType}.}
#'   \item{\code{TERMS} : }{A character vector the same length as the 'setType'
#'      list with the proper names of the sets.}
#'   \item{\code{description} : }{ (OPTIONAL) A character vector the same length
#'      as the 'setType' list with a note on that set (for the \code{.gmt} file
#'      included with this package, this field contains hyperlinks to the
#'      MSigDB description card for that pathway). This field is included when
#'      \code{description = TRUE}.}
#' }
#'
#' @details This function uses \code{R}'s \code{\link{readChar}} function to
#'    improve character input performance over \code{\link{readLines}} (and
#'    far improve input performance over \code{\link{scan}}).
#'
#'    See the Broad Institute's "Data Formats" page for a description of the
#'    Gene Matrix Transposed file format:
#'    \url{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}
#'
#' @export
#'
#' @seealso \code{\link{print.pathwayCollection}}; \code{\link{write_gmt}}
#'
#' @examples
#'   # If you have installed the package:
#'   data_path <- system.file(
#'     "extdata", "c2.cp.v6.0.symbols.gmt",
#'     package = "pathwayPCA", mustWork = TRUE
#'   )
#'   geneset_ls <- read_gmt(data_path, description = TRUE)
#'
#'   # # If you are using the development version from GitHub:
#'   # geneset_ls <- read_gmt(
#'   #   "inst/extdata/c2.cp.v6.0.symbols.gmt",
#'   #   description = TRUE
#'   # )
#'
read_gmt <- function(file, setType = c("pathways", "genes", "regions"),
                     description = FALSE, delim = "\t"){

  # browser()

  # Read the file as a single character vector, split it by line, then split
  #   each line by "tab"
  text_char <- readChar(file, file.info(file)$size, useBytes = TRUE)
  text_vec <- strsplit(text_char, "\n", fixed = TRUE, useBytes = TRUE)[[1]]
  geneset_ls <- strsplit(text_vec, split = delim)

  # Extract the pathway names
  geneset_names <- vapply(geneset_ls, `[[`, 1, FUN.VALUE = character(1))

  # Extract the genes
  genes_ls <- lapply(geneset_ls, function(x){

    x_len <- length(x)
    x[x_len] <- gsub("\r", "", x[x_len])
    x[3:x_len]

  })

  # Create the pathwayCollection output
  out <- list(
    pathways = genes_ls,
    TERMS = geneset_names
  )
  
  setType <- match.arg(setType)
  names(out)[1] <- setType

  if(description){

    geneset_descr <- vapply(geneset_ls, `[[`, 2, FUN.VALUE = character(1))
    out$description <- geneset_descr

  }

  class(out) <- c("pathwayCollection", "list")
  out

}
