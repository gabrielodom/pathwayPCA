% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessClass_pathwayCollection.R
\name{SubsetPathwayCollection}
\alias{SubsetPathwayCollection}
\alias{[[.pathwayCollection}
\title{Subset a \code{pathwayCollection}-class Object by Pathway.}
\usage{
\method{[[}{pathwayCollection}(x, name_char)
}
\arguments{
\item{x}{An object of class \code{pathwayCollection}.}

\item{name_char}{The name of a pathway in the collection or its unique ID.}
}
\value{
A list of the pathway name (\code{Term}), unique ID (\code{pathID}),
   contents (\code{IDs}), description (\code{description}), and number of
   features (\code{Size}).
}
\description{
The subset method for pathways lists as returned by the
   \code{\link{read_gmt}} function.
}
\details{
This function finds the index matching the \code{name_char} argument
   to the \code{TERMS} field of the \code{pathwayCollection}-class Object,
   then subsets the \code{pathways} list, \code{TERMS} vector,
   \code{description} vector, and \code{setsize} vector by this index. If you
   subset a trimmed \code{pathwayCollection} object, and the function errors
   with "Pathway not found.", then the pathway specified has been trimmed
   from the pathway collection.

   Also, this function does not allow for users to overwrite any portion of
   a pathway collection. These objects should rarely, if ever, be changed.
   If you absolutely must change the components of a \code{pathwayCollection}
   object, then create a new one with the code{\link{CreatePathwayCollection}}
   function.
}
\examples{
  data("colon_pathwayCollection")
  colon_pathwayCollection[["KEGG_RETINOL_METABOLISM"]]

}
