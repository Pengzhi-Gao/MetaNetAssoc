test_that("PPI preparation treats reciprocal edges as one edge", {
  ppi <- data.frame(
    left = c("A", "B", "A"),
    right = c("B", "A", "C"),
    evidence = c("db1", "db2", "db1"),
    score = c(800, 700, 600)
  )
  out <- prepare_ppi(
    ppi,
    gene_a = "left",
    gene_b = "right",
    source = "evidence",
    score = "score",
    score_scale = "string"
  )
  expect_equal(nrow(out), 2)
  expect_equal(out$ppi_score[out$gene_a == "A" & out$gene_b == "B"], 0.8)
  expect_match(out$ppi_source[out$gene_a == "A" & out$gene_b == "B"], "db1")
  expect_match(out$ppi_source[out$gene_a == "A" & out$gene_b == "B"], "db2")
})

test_that("wide enrichment member lists are expanded and filtered", {
  x <- data.frame(
    id = c("T1", "T2"),
    description = c("one", "two"),
    members = c("A; B;C", "D;E"),
    fdr = c(0.01, 0.2)
  )
  out <- prepare_enrichment(
    x,
    term_id = "id",
    term_name = "description",
    gene = NULL,
    gene_list = "members",
    p_adjust = "fdr"
  )
  expect_equal(unique(out$term_id), "T1")
  expect_equal(sort(out$gene), c("A", "B", "C"))
})

test_that("invalid scores fail clearly", {
  expect_error(
    prepare_ppi(data.frame(gene_a = "A", gene_b = "B", score = 1200), score = "score"),
    "exceeds 1000"
  )
})

