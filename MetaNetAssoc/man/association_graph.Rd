\name{association_graph}
\alias{association_graph}
\alias{plot_association_network}
\title{Build or Plot a Four-Layer Association Graph}
\usage{
association_graph(x, metabolites = NULL, terms = NULL, min_edge_weight = 0)

plot_association_network(
  x,
  metabolites = NULL,
  terms = NULL,
  min_edge_weight = 0,
  ...
)
}
\arguments{
\item{x}{A \code{MetaNetResult} object.}
\item{metabolites, terms}{Optional IDs or exact names to retain.}
\item{min_edge_weight}{Minimum path weight.}
\item{...}{Additional arguments passed to \code{igraph::plot.igraph}.}
}
\value{An igraph object. The plotting function returns it invisibly.}
\description{Builds the metabolite to MPI gene to term gene to term network.}
