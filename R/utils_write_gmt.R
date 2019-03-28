#' Write a \code{pathwayCollection} Object to a \code{.gmt} File
#'
#' @description Write a \code{pathwayCollection} object as a pathways list file
#'    in Gene Matrix Transposed (\code{.gmt}) format.
#'
#' @param pathwayCollection A \code{pathwayCollection} list of sets. This
#'    list contains the following two or three elements:
#'    \itemize{
#'      \item{'setType' : }{A named list of character vectors. Each vector
#'        contains the names of the individual genes, sites, or CpGs within
#'        that set as a vector of character strings. If you are using genes, 
#'        these genes can be represented by HGNC gene symbols, Entrez IDs,
#'        Ensembl IDs, GO terms, etc.}
#'      \item{\code{TERMS} : }{A character vector the same length as the
#'        'setType' list with the proper names of the sets.}
#'      \item{\code{description} : }{An optional character vector the same
#'        length as the 'setType' list with a note on that set (such as a url
#'        to the description if the set is a pathway). If this element of the
#'        \code{pathwayCollection} is \code{NULL}, then the file will be
#'        written with \code{""} (the empty character string) as its second
#'        field in each line.}
#'    }
#' @param file Either a character string naming a file or a connection open for
#'    writing. File names should end in \code{.gmt} for clarity.
#' @param setType What is the type of the set: pathway set of gene, gene sites
#'    in RNA or DNA, or regions of CpGs. Defaults to \code{''pathway''}.
#'
#' @return NULL. Output written to the file path specified.
#'
#' @details See the Broad Institute's "Data Formats" page for a description of
#'    the Gene Matrix Transposed file format:
#'    \url{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}
#'
#' @export
#'
#' @seealso \code{\link{print.pathwayCollection}}; \code{\link{read_gmt}}
#'
#' @examples
#'   # Toy pathway set
#'   toy_pathwayCollection <- list(
#'     pathways = list(
#'       c("C1orf27", "NR5A1", "BLOC1S4", "C4orf50"),
#'       c("TARS2", "DUSP5", "GPR88"),
#'       c("TRX-CAT3-1", "LINC01333", "LINC01499", "LINC01046", "LINC01149")
#'     ),
#'     TERMS = c("C-or-f_paths", "randomPath2", "randomLINCs"),
#'     description = c("these are", "totally made up", "pathways")
#'   )
#'   class(toy_pathwayCollection) <- c("pathwayCollection", "list")
#'   toy_pathwayCollection
#'
#'   # write_gmt(toy_pathwayCollection, file = "example_pathway.gmt")
#'
write_gmt <- function(pathwayCollection, file,
                      setType = c("pathways", "genes", "regions")){

  ###  Setup  ###
  setType    <- match.arg(setType)
  sets_ls    <- pathwayCollection[setType]
  TERMS_char <- pathwayCollection$TERMS
  desc_char  <- pathwayCollection$description
  nPaths     <- length(sets_ls)

  if(is.null(desc_char)){
    desc_char <- rep("", nPaths)
  }


  ###  Error Checking  ###
  if(nPaths != length(TERMS_char)){
    stop("Number of sets should match number of TERMS.")
  }
  if(nPaths != length(desc_char)){
    stop("Number of sets should match number of set descriptions.")
  }
  sets_ls[which(lengths(sets_ls) == 0)] <- NA_character_


  ###  Write a list in .gmt form  ###
  out_ls <- lapply(seq_len(nPaths), function(i){
    c(TERMS_char[i], desc_char[i], sets_ls[[i]])
  })

  ###  Collapse the list  ###
  out_char <- vapply(out_ls, paste, collapse = "\t", FUN.VALUE = character(1))

  ###  Write the File  ###
  writeLines(out_char, con = file)

}
