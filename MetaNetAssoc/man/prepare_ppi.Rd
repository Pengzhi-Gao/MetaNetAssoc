\name{prepare_ppi}
\alias{prepare_ppi}
\title{Prepare an Undirected Protein-Protein Interaction Table}
\usage{
prepare_ppi(
  x,
  gene_a = "gene_a",
  gene_b = "gene_b",
  source = NULL,
  score = NULL,
  score_scale = c("auto", "unit", "string")
)
}
\arguments{
\item{x}{A data frame containing PPI edges.}
\item{gene_a, gene_b, source, score}{Column names.}
\item{score_scale}{Automatic detection, an existing 0-1 scale, or STRING's
0-1000 scale.}
}
\value{A standardized, deduplicated data frame.}
\description{Reciprocal and repeated PPI edges are collapsed to avoid double counting.}

