\name{load_mpi_reference}
\alias{load_mpi_reference}
\title{Load a packaged human or mouse MPI reference network}
\usage{
load_mpi_reference(species = c("human", "mouse"), include_metadata = FALSE)
}
\arguments{
  \item{species}{Either \code{"human"} or \code{"mouse"}.}
  \item{include_metadata}{Whether to return the build metadata together with the network.}
}
\value{A data frame, or a list with \code{network} and \code{metadata}.}
\description{Loads the bundled MPI references. The human network integrates
project-provided KEGG, Reactome, Human-GEM, and BRENDA evidence. The mouse
network is an Ensembl Compara orthology projection of the human evidence. Both
use KEGG compound IDs, not HMDB IDs.}
\examples{
human_mpi <- load_mpi_reference("human")
head(human_mpi)
}
