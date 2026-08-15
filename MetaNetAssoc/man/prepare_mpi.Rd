\name{prepare_mpi}
\alias{prepare_mpi}
\title{Prepare a Metabolite-Protein Interaction Table}
\usage{
prepare_mpi(
  x,
  metabolite_id = "metabolite_id",
  metabolite_name = "metabolite_name",
  gene = "relation_gene",
  source = NULL,
  weight = NULL,
  weight_scale = c("auto", "unit", "string")
)
}
\arguments{
\item{x}{A data frame containing metabolite-gene associations.}
\item{metabolite_id, metabolite_name, gene, source, weight}{Column names. At least
one metabolite ID or name column must be supplied.}
\item{weight_scale}{Automatic detection, an existing 0-1 scale, or a 0-1000 scale.}
}
\value{A standardized data frame.}
\description{Validates and standardizes metabolite-gene evidence.}

