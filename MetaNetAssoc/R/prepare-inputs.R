#' Prepare a metabolite-protein interaction table
#'
#' @param x A data frame containing metabolite-gene associations.
#' @param metabolite_id,metabolite_name,gene,source,weight Column names.
#'   At least one of `metabolite_id` and `metabolite_name` must be supplied.
#' @param weight_scale Either automatic detection, an existing 0-1 scale, or
#'   STRING-like 0-1000 values.
#' @return A standardized data frame.
#' @export
prepare_mpi <- function(
    x,
    metabolite_id = "metabolite_id",
    metabolite_name = "metabolite_name",
    gene = "relation_gene",
    source = NULL,
    weight = NULL,
    weight_scale = c("auto", "unit", "string")) {
  .assert_data_frame(x, "x")
  weight_scale <- match.arg(weight_scale)
  if (is.null(metabolite_id) && is.null(metabolite_name)) {
    stop("Supply at least one metabolite ID or name column.", call. = FALSE)
  }
  .assert_columns(x, gene, "x")

  id <- .clean_text(.column_or(x, metabolite_id, NA_character_))
  name <- .clean_text(.column_or(x, metabolite_name, NA_character_))
  relation_gene <- .clean_gene(x[[gene]])
  mpi_source <- .clean_text(.column_or(x, source, "unspecified"))
  mpi_source[is.na(mpi_source)] <- "unspecified"
  mpi_weight <- if (is.null(weight)) rep(1, nrow(x)) else x[[weight]]
  mpi_weight <- .normalise_unit_score(mpi_weight, "MPI weight", weight_scale)

  keep <- (!is.na(id) | !is.na(name)) & !is.na(relation_gene)
  out <- data.frame(
    metabolite_id = id[keep],
    metabolite_name = name[keep],
    relation_gene = relation_gene[keep],
    mpi_source = mpi_source[keep],
    mpi_weight = mpi_weight[keep],
    stringsAsFactors = FALSE
  )
  unique(out)
}

#' Prepare an undirected protein-protein interaction table
#'
#' Reciprocal and repeated edges are collapsed so that a PPI is never counted
#' twice. For repeated edges the maximum confidence is retained and source names
#' are concatenated.
#'
#' @param x A data frame containing PPI edges.
#' @param gene_a,gene_b,source,score Column names.
#' @param score_scale Either automatic detection, an existing 0-1 scale, or
#'   STRING's 0-1000 range.
#' @return A standardized, deduplicated data frame.
#' @export
prepare_ppi <- function(
    x,
    gene_a = "gene_a",
    gene_b = "gene_b",
    source = NULL,
    score = NULL,
    score_scale = c("auto", "unit", "string")) {
  .assert_data_frame(x, "x")
  .assert_columns(x, c(gene_a, gene_b), "x")
  score_scale <- match.arg(score_scale)

  a <- .clean_gene(x[[gene_a]])
  b <- .clean_gene(x[[gene_b]])
  ppi_source <- .clean_text(.column_or(x, source, "unspecified"))
  ppi_source[is.na(ppi_source)] <- "unspecified"
  ppi_score <- if (is.null(score)) rep(1, nrow(x)) else x[[score]]
  ppi_score <- .normalise_unit_score(ppi_score, "PPI score", score_scale)

  keep <- !is.na(a) & !is.na(b) & a != b
  a <- a[keep]
  b <- b[keep]
  ppi_source <- ppi_source[keep]
  ppi_score <- ppi_score[keep]
  first <- ifelse(a <= b, a, b)
  second <- ifelse(a <= b, b, a)
  key <- paste(first, second, sep = "\r")
  groups <- split(seq_along(key), key)

  if (!length(groups)) {
    return(data.frame(
      gene_a = character(), gene_b = character(),
      ppi_source = character(), ppi_score = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, lapply(groups, function(i) {
    data.frame(
      gene_a = first[i[1]],
      gene_b = second[i[1]],
      ppi_source = .collapse_unique(ppi_source[i]),
      ppi_score = max(ppi_score[i]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' Prepare enrichment results as one term-gene pair per row
#'
#' @param x A data frame containing enrichment results.
#' @param term_id,term_name,gene,gene_list,group,p_value,p_adjust Column names.
#'   Use `gene` for long data or `gene_list` for a delimited member column.
#' @param gene_sep Regular expression used to split `gene_list`.
#' @param fdr_cutoff,pvalue_cutoff Optional significance thresholds. FDR is used
#'   when an adjusted p-value column is available; otherwise p-value is used.
#' @return A standardized long enrichment table.
#' @export
prepare_enrichment <- function(
    x,
    term_id = "term_id",
    term_name = "term_name",
    gene = "gene",
    gene_list = NULL,
    group = NULL,
    p_value = NULL,
    p_adjust = NULL,
    gene_sep = "[;,/|[:space:]]+",
    fdr_cutoff = 0.05,
    pvalue_cutoff = 0.05) {
  .assert_data_frame(x, "x")
  .assert_columns(x, c(term_id, term_name), "x")
  if (is.null(gene) && is.null(gene_list)) {
    stop("Supply either `gene` or `gene_list`.", call. = FALSE)
  }
  if (!is.null(gene)) .assert_columns(x, gene, "x")
  if (!is.null(gene_list)) .assert_columns(x, gene_list, "x")

  n <- nrow(x)
  term_id_value <- .clean_text(x[[term_id]])
  term_name_value <- .clean_text(x[[term_name]])
  group_value <- .clean_text(.column_or(x, group, "Member", n))
  group_value[is.na(group_value)] <- "Member"
  p_value_value <- suppressWarnings(as.numeric(.column_or(x, p_value, NA_real_, n)))
  p_adjust_value <- suppressWarnings(as.numeric(.column_or(x, p_adjust, NA_real_, n)))

  keep <- !is.na(term_id_value) & !is.na(term_name_value)
  if (!is.null(p_adjust) && !is.null(fdr_cutoff)) {
    keep <- keep & !is.na(p_adjust_value) & p_adjust_value <= fdr_cutoff
  } else if (!is.null(p_value) && !is.null(pvalue_cutoff)) {
    keep <- keep & !is.na(p_value_value) & p_value_value <= pvalue_cutoff
  }

  x <- x[keep, , drop = FALSE]
  term_id_value <- term_id_value[keep]
  term_name_value <- term_name_value[keep]
  group_value <- group_value[keep]
  p_value_value <- p_value_value[keep]
  p_adjust_value <- p_adjust_value[keep]

  genes <- if (!is.null(gene)) {
    lapply(.clean_gene(x[[gene]]), function(z) z[!is.na(z)])
  } else {
    .split_rows(x[[gene_list]], gene_sep)
  }
  lengths <- lengths(genes)
  if (!sum(lengths)) {
    stop("No valid term genes remained after preparation.", call. = FALSE)
  }

  out <- data.frame(
    term_id = rep(term_id_value, lengths),
    term_name = rep(term_name_value, lengths),
    gene = unlist(genes, use.names = FALSE),
    group = rep(group_value, lengths),
    p_value = rep(p_value_value, lengths),
    p_adjust = rep(p_adjust_value, lengths),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene), , drop = FALSE]
  out <- out[!duplicated(out[c("term_id", "gene")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

