test_that("synthetic example runs end to end", {
  example <- load_example_data()
  fit <- run_metanet(
    example$metabolites$metabolite_id,
    mpi = example$mpi,
    ppi = example$ppi,
    enrichment = example$enrichment,
    focus_genes = example$focus_genes
  )

  expect_s4_class(fit, "MetaNetResult")
  expect_equal(nrow(result_scores(fit)), 10 * 7)
  expect_gt(nrow(result_edges(fit)), 0)
  expect_true(all(result_scores(fit)$normalized_score >= 0))
  expect_true(all(result_scores(fit)$normalized_score <= 1))
  expect_true(all(c("ppi", "gene_overlap") %in% unique(result_edges(fit)$path_type)))
})

test_that("dominant tutorial signals rank first for representative terms", {
  example <- load_example_data()
  fit <- run_metanet(
    example$metabolites$metabolite_id,
    mpi = example$mpi,
    ppi = example$ppi,
    enrichment = example$enrichment
  )
  scores <- result_scores(fit)
  top_name <- function(term) {
    rows <- scores[scores$term_id == term, ]
    rows$metabolite_name[which.max(rows$normalized_score)]
  }
  expect_equal(top_name("GO:SIM0003"), "Glutathione")
  expect_true(top_name("GO:SIM0004") %in% c("L-Carnitine", "Acetylcarnitine", "Pyruvate"))
  expect_equal(top_name("GO:SIM0006"), "Spermidine")
})

test_that("each tutorial metabolite recovers its designed primary term", {
  example <- load_example_data()
  fit <- run_metanet(
    example$metabolites$metabolite_id,
    mpi = example$mpi,
    ppi = example$ppi,
    enrichment = example$enrichment
  )
  scores <- result_scores(fit)
  top <- do.call(rbind, lapply(split(scores, scores$metabolite_id), function(rows) {
    rows[which.max(rows$normalized_score), c("metabolite_id", "term_id")]
  }))
  expected <- example$metabolites[c("metabolite_id", "expected_primary_term")]
  observed <- merge(expected, top, by = "metabolite_id", sort = FALSE)
  expect_equal(observed$term_id, observed$expected_primary_term)
})

test_that("unmatched metabolites are reported without losing matched inputs", {
  example <- load_example_data()
  fit <- run_metanet(
    c("SIMMET0001", "not-present"),
    mpi = example$mpi,
    ppi = example$ppi,
    enrichment = example$enrichment
  )
  expect_equal(result_unmatched(fit)$metabolites, "not-present")
  expect_equal(length(unique(result_scores(fit)$metabolite_id)), 1)
})

test_that("hub modes retain bounded scores", {
  example <- load_example_data()
  fits <- lapply(c("none", "favor", "penalize"), function(mode) {
    run_metanet(
      example$metabolites$metabolite_id,
      mpi = example$mpi,
      ppi = example$ppi,
      enrichment = example$enrichment,
      hub_mode = mode
    )
  })
  expect_true(all(vapply(fits, function(x) max(result_scores(x)$normalized_score) <= 1, logical(1))))
})

test_that("PPI source filtering is applied before reciprocal-edge collapse", {
  example <- load_example_data()
  fit <- run_metanet(
    "SIMMET0001",
    mpi = example$mpi,
    ppi = example$ppi,
    enrichment = example$enrichment,
    ppi_sources = "STRING_mock"
  )
  ppi_paths <- result_edges(fit)
  ppi_paths <- ppi_paths[ppi_paths$path_type == "ppi", ]
  expect_true(all(ppi_paths$ppi_source == "STRING_mock"))
})
