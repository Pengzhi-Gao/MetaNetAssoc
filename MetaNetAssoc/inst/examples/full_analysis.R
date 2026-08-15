library(MetaNetAssoc)

example <- load_example_data()

fit <- run_metanet(
  metabolites = example$metabolites$metabolite_id,
  focus_genes = example$focus_genes,
  mpi = example$mpi,
  ppi = example$ppi,
  enrichment = example$enrichment,
  include_overlap = TRUE,
  hub_mode = "none"
)

print(fit)
print(utils::head(result_scores(fit), 10))
print(result_qc(fit))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  print(plot_term_ranking(fit, "Response to oxidative stress"))
  print(plot_score_heatmap(fit))
}

