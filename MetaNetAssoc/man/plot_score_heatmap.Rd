\name{plot_score_heatmap}
\alias{plot_score_heatmap}
\title{Plot a Metabolite-Term Score Heatmap}
\usage{
plot_score_heatmap(
  x,
  metric = "normalized_score",
  top_n_metabolites = NULL,
  top_n_terms = NULL
)
}
\arguments{
\item{x}{A \code{MetaNetResult} object.}
\item{metric}{Numeric score column.}
\item{top_n_metabolites, top_n_terms}{Optional display limits.}
}
\value{A ggplot object.}
\description{Draws the score matrix as a heatmap.}

