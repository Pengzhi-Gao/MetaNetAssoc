\name{run_metanet}
\alias{run_metanet}
\title{Associate Metabolites with Biological Terms}
\usage{
run_metanet(
  metabolites = NULL,
  mpi = NULL,
  ppi = NULL,
  enrichment = NULL,
  focus_genes = NULL,
  mpi_sources = NULL,
  ppi_sources = NULL,
  include_overlap = TRUE,
  hub_mode = c("none", "favor", "penalize"),
  term_similarity_cutoff = NULL,
  keep_zero = TRUE,
  deg = NULL,
  species = c("human", "mouse"),
  mpi_database = NULL,
  string_data = NULL,
  string_score_threshold = 400
)
}
\arguments{
\item{metabolites}{Candidate metabolite IDs and/or names.}
\item{mpi}{Standardized MPI data from \code{prepare_mpi()}.}
\item{ppi}{Standardized PPI data from \code{prepare_ppi()}.}
\item{enrichment}{Standardized term-gene data.}
\item{focus_genes}{Optional original gene set used for QC reporting.}
\item{mpi_sources, ppi_sources}{Optional source names to retain.}
\item{include_overlap}{Whether direct MPI-gene membership in a term is a path.}
\item{hub_mode}{No hub adjustment, favor term-gene hubs, or penalize hubs.}
\item{term_similarity_cutoff}{Optional Jaccard threshold for term reduction.}
\item{keep_zero}{Whether to retain pairs with no supporting path.}
\item{deg}{Optional differential genes used to build missing input networks.}
\item{species}{\code{"human"} or \code{"mouse"}.}
\item{mpi_database}{HMDB--gene evidence table for automatic MPI construction.}
\item{string_data}{Optional offline STRING edge table.}
\item{string_score_threshold}{Minimum STRING combined score.}
}
\value{A \code{MetaNetResult} object.}
\description{
Finds direct links between each metabolite's MPI genes and each term's genes,
retains supporting paths, and calculates raw and normalized association scores.
}
\examples{
example <- load_example_data()
fit <- run_metanet(
  example$metabolites$metabolite_id,
  mpi = example$mpi,
  ppi = example$ppi,
  enrichment = example$enrichment
)
head(result_scores(fit))
}
