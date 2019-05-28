#' Wikipathways Homosapiens Gene Symbols
#'
#' @description A \code{pathwayCollection} object containing the homosapiens
#'    pathways list from Wikipathways (\url{https://www.wikipathways.org/}).
#'
#' @details This \code{pathwayCollection} was sent to us from Dr. Alexander Pico
#'    at the Gladstone Institute
#'    (\url{https://gladstone.org/our-science/people/alexander-pico}).
#'    
#'    This pathway collection was translated from EntrezIDs to HGNC Symbols with
#'    the script \code{convert_EntrezID_to_HGNC_Ensembl.R} in \code{scripts}.
#'
#' @format A \code{pathwayCollection} list of three elements:
#' \itemize{
#'   \item{\code{pathways} : }{A named list of 457 character vectors. Each
#'      vector contains the Gene Symbols of the individual genes within that
#'      pathway as a vector of character strings. The names are the shorthand
#'      pathway names.}
#'   \item{\code{TERMS} : }{A character vector of length 457 containing the
#'      shorthand names of the gene pathways.}
#'   \item{\code{description} : }{A character vector of length 457 containing
#'      the full names of the gene pathways.}
#' }
#'
#' @source Dr. Alexander Pico, Wikipathways
"wikipwsHS_Symbol_pathwayCollection"
