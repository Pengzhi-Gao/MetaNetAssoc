test_that("term reduction keeps the more significant representative", {
  enrichment <- data.frame(
    term_id = rep(c("T1", "T2", "T3"), each = 3),
    term_name = rep(c("first", "second", "third"), each = 3),
    gene = c("A", "B", "C", "A", "B", "C", "X", "Y", "Z"),
    group = "G",
    p_value = rep(c(0.001, 0.01, 0.02), each = 3),
    p_adjust = rep(c(0.003, 0.02, 0.03), each = 3)
  )
  reduced <- reduce_terms(enrichment, jaccard_cutoff = 0.8)
  expect_equal(sort(unique(reduced$enrichment$term_id)), c("T1", "T3"))
  t2 <- reduced$mapping[reduced$mapping$term_id == "T2", ]
  expect_equal(t2$representative_term_id, "T1")
  expect_false(t2$retained)
})

