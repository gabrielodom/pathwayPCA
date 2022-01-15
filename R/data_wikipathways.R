#' Wikipathways Homosapiens EntrezIDs
#'
#' @description A \code{pathwayCollection} object containing the homosapiens
#'    pathways list from Wikipathways (\url{https://www.wikipathways.org/}).
#'
#' @details This \code{pathwayCollection} was sent to us from Dr. Alexander Pico
#'    at the Gladstone Institute
#'    (\url{https://gladstone.org/our-science/people/alexander-pico}).
#'
#' @format A \code{pathwayCollection} list of three elements:
#' \itemize{
#'   \item{\code{pathways} : }{A named list of 443 character vectors. Each
#'      vector contains the Entrez Gene IDs of the individual genes within that
#'      pathway as a vector of character strings. The names are the shorthand
#'      pathway names.}
#'   \item{\code{TERMS} : }{A character vector of length 443 containing the
#'      shorthand names of the gene pathways.}
#'   \item{\code{description} : }{A character vector of length 443 containing
#'      the full names of the gene pathways.}
#' }
#'
#' @source Dr. Alexander Pico, Wikipathways
#' @usage data(wikipwsHS_Entrez_pathwayCollection)
"wikipwsHS_Entrez_pathwayCollection"
