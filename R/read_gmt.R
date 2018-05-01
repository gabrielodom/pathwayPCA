#' Read a \code{.gmt} file in as a \code{pathwaySet} object
#'
#' @description Read a pathways list file in Gene Matrix Transposed
#'    (\code{.gmt}) format, with special performance consideration for large
#'    files. Present this object as a \code{pathwaySet} object.
#'
#' @param file A path to a file or a connection. This file must be a \code{.gmt}
#'    file, otherwise input will likely be nonsense. See the "Details" section
#'    for more information.
#' @param description Should the "description" field (the second field in the
#'    \code{.gmt} file on each line) be included in the output? Defaults to
#'    \code{FALSE}.
#' @param delim The \code{.gmt} delimiter. As proper \code{.gmt} files are tab
#'    delimited, this defaults to \code{"\\t"}.
#'
#' @return A \code{pathwaySet} list of pathways. This list has three elements:
#' \itemize{
#'   \item{\code{pathways} : }{A named list of character vectors. Each vector
#'      contains the names of the individual genes within that pathway as a
#'      vector of character strings.}
#'   \item{\code{TERMS} : }{A character vector the same length as the
#'      \code{pathways} list with the proper names of the pathways.}
#'   \item{\code{description} : }{ (OPTIONAL) A character vector the same length
#'      as the \code{pathways} list with a note on that pathway (for the
#'      \code{.gmt} file included with this package, this field contains
#'      hyperlinks to the MSigDB description card for that pathway). This field
#'      is included when \code{description = TRUE}.}
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
#' @seealso \code{\link{print.pathwaySet}}
#'
#' @examples
#' \dontrun{
#'   # If you have installed the package:
#'   data_path <- system.file("extdata", "c2.cp.v6.0.symbols.gmt",
#'                             package = "pathwayPCA", mustWork = TRUE)
#'   geneset_ls <- read_gmt(data_path)
#'
#'   # If you are using the development version from GitHub:
#'   geneset_ls <- read_gmt("inst/extdata/c2.cp.v6.0.symbols.gmt")
#' }
#'
read_gmt <- function(file, description = FALSE, delim = "\t"){

  # Read the file as a single character vector, split it by line, then split
  #   each line by "tab"
  text_char <- readChar(file, file.info(file)$size, useBytes = TRUE)
  text_vec <- strsplit(text_char, "\n", fixed = TRUE, useBytes = TRUE)[[1]]
  geneset_ls <- strsplit(text_vec, split = delim)

  # Extract the pathway names
  geneset_names <- sapply(geneset_ls, `[[`, 1)

  # Extract the genes
  genes_ls <- lapply(geneset_ls, function(x){

    x_len <- length(x)
    x[3:x_len]

  })

  # Create the pathwaySet output
  out <- list(pathways = genes_ls,
              TERMS = geneset_names)

  if(description){

    geneset_descr <- sapply(geneset_ls, `[[`, 2)
    out$description <- geneset_descr

  }

  class(out) <- c("pathwaySet", "list")
  out

}
