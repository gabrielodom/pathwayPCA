#' Write a \code{pathwaySet} Object to a \code{.gmt} File
#'
#' @description Write a \code{pathwaySet} object as a pathways list file in
#'    Gene Matrix Transposed (\code{.gmt}) format.
#'
#' @param pathwaySet A \code{pathwaySet} list of pathways. This list contains
#'    the following three elements:
#'    \itemize{
#'      \item{\code{pathways} : }{A named list of character vectors. Each vector
#'        contains the names of the individual genes within that pathway as a
#'        vector of character strings. Genes can be represented by HGNC gene
#'        symbols, Entrez IDs, Ensembl IDs, GO terms, etc. }
#'      \item{\code{TERMS} : }{A character vector the same length as the
#'        \code{pathways} list with the proper names of the pathways.}
#'      \item{\code{description} : }{A character vector the same length as the
#'        \code{pathways} list with a note on that pathway (such as a url to the
#'        description of the pathway). If this element of \code{pathwaySet} is
#'        \code{NULL}, then the file will be written with \code{""} (the empty
#'        character string) as its second field in each line.}
#'    }
#' @param file Either a character string naming a file or a connection open for
#'    writing. File names should end in \code{.gmt} for clarity
#'
#' @details See the Broad Institute's "Data Formats" page for a description of
#'    the Gene Matrix Transposed file format:
#'    \url{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}
#'
#' @export
#'
#' @seealso \code{\link{print.pathwaySet}}; \code{\link{read_gmt}}
#'
#' @examples
#' \dontrun{
#'   # Toy pathway set
#'   toy_pathwaySet <- list(
#'     pathways = list(
#'       c("C1orf27", "NR5A1", "BLOC1S4", "C4orf50"),
#'       c("TARS2", "DUSP5", "GPR88"),
#'       c("TRX-CAT3-1", "LINC01333", "LINC01499", "LINC01046", "LINC01149")
#'     ),
#'     TERMS = c("C-or-f_paths", "randomPath2", "randomLINCs"),
#'     description = c("these are", "totally made up", "pathways")
#'   )
#'   class(toy_pathwaySet) <- c("pathwaySet", "list")
#'   print(toy_pathwaySet)
#'
#'   write_gmt(toy_pathwaySet, file = "example_pathway.gmt")
#' }
#'
write_gmt <- function(pathwaySet, file){

  ###  Setup  ###
  pathways_ls <- pathwaySet$pathways
  TERMS_char  <- pathwaySet$TERMS
  desc_char   <- pathwaySet$description
  nPaths      <- length(pathways_ls)

  if(is.null(desc_char)){
    desc_char <- rep("", nPaths)
  }


  ###  Error Checking  ###
  if(nPaths != length(TERMS_char)){
    stop("Number of pathways should match number of TERMS.")
  }
  if(nPaths != length(desc_char)){
    stop("Number of pathways should match number of pathway descriptions.")
  }
  pathways_ls[which(lengths(pathways_ls) == 0)] <- NA_character_


  ###  Write a list in .gmt form  ###
  out_ls <- lapply(1:nPaths, function(i){
    c(TERMS_char[i], desc_char[i], pathways_ls[[i]])
  })

  ###  Collapse the list  ###
  out_char <- sapply(out_ls, paste, collapse = "\t")

  ###  Write the File  ###
  writeLines(out_char, con = file)

}
