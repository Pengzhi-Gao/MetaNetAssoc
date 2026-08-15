\name{reduce_terms}
\alias{reduce_terms}
\title{Reduce Redundant Terms}
\usage{
reduce_terms(enrichment, jaccard_cutoff = 0.7, within_group = TRUE)
}
\arguments{
\item{enrichment}{A standardized long enrichment table.}
\item{jaccard_cutoff}{Member-gene Jaccard threshold between zero and one.}
\item{within_group}{Whether to compare terms only within the same group.}
}
\value{A list containing reduced \code{enrichment} and a \code{mapping} table.}
\description{Greedily keeps the strongest term among highly similar member-gene sets.}

