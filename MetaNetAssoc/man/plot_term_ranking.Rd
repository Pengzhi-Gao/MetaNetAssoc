\name{plot_term_ranking}
\alias{plot_term_ranking}
\title{Plot the Metabolite Ranking for One Term}
\usage{
plot_term_ranking(x, term, metric = "normalized_score", top_n = 20)
}
\arguments{
\item{x}{A \code{MetaNetResult} object.}
\item{term}{A term ID or exact name.}
\item{metric}{Numeric score column.}
\item{top_n}{Maximum metabolites to show.}
}
\value{A ggplot object.}
\description{Draws a horizontal metabolite ranking for a selected term.}

