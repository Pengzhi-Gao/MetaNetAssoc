\name{MetaNetResult-class}
\alias{MetaNetResult-class}
\alias{show,MetaNetResult-method}
\title{MetaNetAssoc Result Class}
\description{An S4 container returned by \code{run_metanet()}.}
\section{Slots}{
\describe{
  \item{\code{call}}{The matched analysis call.}
  \item{\code{parameters}}{Analysis settings.}
  \item{\code{scores}}{Metabolite-term score table.}
  \item{\code{edges}}{Auditable metabolite-gene-gene-term paths.}
  \item{\code{term_table}}{Terms retained for analysis.}
  \item{\code{qc}}{Quality-control summary.}
  \item{\code{unmatched}}{Unmatched metabolites and genes.}
  \item{\code{term_reduction}}{Optional redundant-term mapping.}
}}

