\name{prepare_enrichment}
\alias{prepare_enrichment}
\title{Prepare Enrichment Results}
\usage{
prepare_enrichment(
  x,
  term_id = "term_id",
  term_name = "term_name",
  gene = "gene",
  gene_list = NULL,
  group = NULL,
  p_value = NULL,
  p_adjust = NULL,
  gene_sep = "[;,/|[:space:]]+",
  fdr_cutoff = 0.05,
  pvalue_cutoff = 0.05
)
}
\arguments{
\item{x}{A data frame containing enrichment results.}
\item{term_id, term_name, gene, gene_list, group, p_value, p_adjust}{Column names.}
\item{gene_sep}{Regular expression used to split a member-gene list.}
\item{fdr_cutoff, pvalue_cutoff}{Optional significance thresholds.}
}
\value{A standardized long table with one term-gene pair per row.}
\description{Converts long or delimited enrichment results into the internal format.}

