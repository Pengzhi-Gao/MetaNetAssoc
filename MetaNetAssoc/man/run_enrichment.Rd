\name{run_enrichment}
\alias{run_enrichment}
\title{Run GO Enrichment with clusterProfiler}
\usage{
run_enrichment(
  genes,
  organism = c("human", "mouse"),
  ontology = "BP",
  background = NULL,
  p_adjust_cutoff = 0.05,
  min_gene_set_size = 10,
  max_gene_set_size = 500
)
}
\arguments{
\item{genes}{Gene symbols.}
\item{organism}{Human or mouse.}
\item{ontology}{GO ontology passed to \code{clusterProfiler::enrichGO}.}
\item{background}{Optional background gene symbols.}
\item{p_adjust_cutoff}{Adjusted p-value cutoff.}
\item{min_gene_set_size, max_gene_set_size}{Term size limits.}
}
\value{A standardized long enrichment table.}
\description{An optional convenience wrapper around clusterProfiler.}

