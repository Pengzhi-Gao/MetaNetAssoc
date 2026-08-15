\name{result_accessors}
\alias{result_accessors}
\alias{result_scores}
\alias{result_edges}
\alias{result_qc}
\alias{result_unmatched}
\alias{result_term_reduction}
\title{Extract Tables from a MetaNetAssoc Result}
\usage{
result_scores(x)

result_edges(x)

result_qc(x)

result_unmatched(x)

result_term_reduction(x)
}
\arguments{\item{x}{A \code{MetaNetResult} object.}}
\value{A data frame, except \code{result_unmatched()} which returns a list.}
\description{Access result components without direct S4 slot manipulation.}

