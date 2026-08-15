test_that("ggplot helpers return plots", {
  skip_if_not_installed("ggplot2")
  example <- load_example_data()
  fit <- run_metanet(
    example$metabolites$metabolite_id,
    mpi = example$mpi,
    ppi = example$ppi,
    enrichment = example$enrichment
  )
  expect_s3_class(plot_term_ranking(fit, "GO:SIM0003"), "ggplot")
  expect_s3_class(plot_score_heatmap(fit), "ggplot")
})
