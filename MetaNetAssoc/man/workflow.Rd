\name{workflow}
\alias{build_ppi_network}
\alias{build_mpi_network}
\alias{run_deg_enrichment}
\alias{export_enrichment}
\title{DEG and HMDB workflow helpers}
\usage{
build_ppi_network(deg, species = c("human", "mouse"), gene_column = "gene",
  score_threshold = 400, string_data = NULL, string_version = "12-0")
build_mpi_network(deg, metabolites, species = c("human", "mouse"),
  mpi_database = NULL, gene_column = "gene", hmdb_column = "hmdb_id")
run_deg_enrichment(deg, species = c("human", "mouse"), gene_column = "gene",
  export_file = NULL, ...)
export_enrichment(enrichment, file)
}
\arguments{
  \item{deg}{A gene vector or data frame.}
  \item{metabolites}{HMDB identifiers or a data frame with an HMDB-ID column.}
  \item{species}{\code{"human"} or \code{"mouse"}.}
  \item{gene_column,hmdb_column}{Column names for data-frame inputs.}
  \item{score_threshold}{Minimum STRING combined score, 0--1000.}
  \item{string_data}{Optional downloaded STRING edge table for offline use.}
  \item{string_version}{STRING API version, e.g. \code{"12-0"} or \code{"12.0"}.}
  \item{mpi_database}{An evidence table mapping HMDB IDs to gene symbols.}
  \item{export_file,file}{CSV output path.}
  \item{enrichment}{A result table or \code{run_deg_enrichment()} result.}
  \item{...}{Passed to \code{run_enrichment}.}
}
\value{PPI and MPI builders return standardized edge tables. Enrichment returns
a list containing raw and normalized tables. \code{export_enrichment} invisibly
returns the output path.}
\description{Helpers implementing the DEG/HMDB analysis workflow. STRING is
queried only for PPI. MPI evidence must be supplied as a traceable reference
table; the bundled reference is synthetic and is suitable only for examples.}
\examples{
dat <- load_example_data()
mpi <- build_mpi_network(c("CHKA", "SLC44A1"), "HMDB0000122",
  mpi_database = dat$hmdb_mpi_reference)
}
