# MetaNetAssoc 0.2.0: fully offline DEG + HMDB workflow.
# The bundled MPI evidence is synthetic; replace it with a documented
# HMDB--gene evidence table for biological analyses.

library(MetaNetAssoc)

dat <- load_example_data()
deg <- c("CHKA", "SLC44A1")
hmdb <- "HMDB0000122"

# 1. Use a saved STRING-like table. Scores here are on STRING's 0--1000 scale.
string_offline <- data.frame(
  gene_a = dat$ppi$gene_a,
  gene_b = dat$ppi$gene_b,
  combined_score = dat$ppi$ppi_score * 1000,
  stringsAsFactors = FALSE
)
ppi <- build_ppi_network(
  deg,
  species = "human",
  string_data = string_offline,
  score_threshold = 400
)

# 2. Construct MPI from DEG, HMDB IDs, and an evidence table.
mpi <- build_mpi_network(
  deg,
  hmdb,
  species = "human",
  mpi_database = dat$hmdb_mpi_reference
)

# 3. For a real analysis, replace this prepared demonstration table with:
# enr <- run_deg_enrichment(deg, species = "human", export_file = "go.csv")
# enrichment <- enr$enrichment
enrichment <- dat$enrichment

# 4. Integrate and score.
fit <- run_metanet(
  metabolites = hmdb,
  mpi = mpi,
  ppi = ppi,
  enrichment = enrichment,
  focus_genes = deg
)

print(fit)
print(result_scores(fit))

# 5. Optional plotting (requires ggplot2; network plot also requires igraph).
# plot_term_ranking(fit)
# plot_score_heatmap(fit)
# plot_association_network(fit)
