% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreatePathwayCollection.R
\name{CreatePathwayCollection}
\alias{CreatePathwayCollection}
\title{Manually Create a \code{pathwayCollection}-class Object.}
\usage{
CreatePathwayCollection(
  sets_ls,
  TERMS,
  setType = c("pathways", "genes", "regions"),
  ...
)
}
\arguments{
\item{sets_ls}{A named list of character vectors. Each vector should contain
the names of the individual genes, proteins, sits, or CpGs within that
set as a vector of character strings. If you create this pathway
collection to integrate with data of class \code{Omics}*, the names
contained in these vectors should have non-empty overlap with the feature
names of the assay data frame that will be paired with this list in the
subsequent analysis.}

\item{TERMS}{A character vector the same length as the \code{sets_ls} list
with the proper names of the sets.}

\item{setType}{What is the type of the set: pathway set of gene, gene sites
in RNA or DNA, or regions of CpGs. Defaults to \code{''pathway''}.}

\item{...}{Additional vectors or data components related to the
\code{sets_ls} list. These values should be passed as a name-value pair.
See "Details" for more information.}
}
\value{
A list object with class \code{pathwayCollection}.
}
\description{
Manually create a \code{pathwayCollection} list similar to the
   output of the \code{\link{read_gmt}} function.
}
\details{
This function checks the set list and set term inputs and then
   creates a \code{pathwayCollection} object from them. Pass additional
   list elements (such as the \code{description} of each set) using the
   form \code{tag = value} through the \code{...} argument (as in the
   \code{\link{list}} function). Because some functions in the
   \code{pathwayPCA} package add and edit elements of \code{pathwayCollection}
   objects, please do not create \code{pathwayCollection} list items named
   \code{setsize} or \code{n_tested}.
}
\examples{
  data("colon_pathwayCollection")

  CreatePathwayCollection(
    sets_ls = colon_pathwayCollection$pathways,
    TERMS = colon_pathwayCollection$TERMS
  )

}
\seealso{
\code{\link{read_gmt}}
}
