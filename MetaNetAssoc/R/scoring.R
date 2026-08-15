.aggregate_mpi <- function(mpi) {
  met_key <- ifelse(
    !is.na(mpi$metabolite_id),
    paste0("id:", toupper(mpi$metabolite_id)),
    paste0("name:", tolower(mpi$metabolite_name))
  )
  key <- paste(met_key, mpi$relation_gene, sep = "\r")
  groups <- split(seq_len(nrow(mpi)), key)
  out <- do.call(rbind, lapply(groups, function(i) {
    weights <- mpi$mpi_weight[i]
    data.frame(
      metabolite_key = met_key[i[1]],
      metabolite_id = mpi$metabolite_id[i[1]],
      metabolite_name = mpi$metabolite_name[i[1]],
      relation_gene = mpi$relation_gene[i[1]],
      mpi_source = .collapse_unique(mpi$mpi_source[i]),
      mpi_weight = 1 - prod(1 - weights),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

.empty_path_table <- function() {
  data.frame(
    metabolite_key = character(),
    metabolite_id = character(),
    metabolite_name = character(),
    relation_gene = character(),
    mpi_source = character(),
    mpi_weight = numeric(),
    enrich_gene = character(),
    ppi_source = character(),
    ppi_score = numeric(),
    term_id = character(),
    term_name = character(),
    group = character(),
    p_value = numeric(),
    p_adjust = numeric(),
    path_type = character(),
    term_gene_degree = integer(),
    hub_factor = numeric(),
    edge_weight = numeric(),
    stringsAsFactors = FALSE
  )
}

.build_paths <- function(mpi, ppi, enrichment, include_overlap, hub_mode) {
  common_columns <- c(
    "metabolite_key", "metabolite_id", "metabolite_name", "relation_gene", "mpi_source",
    "mpi_weight", "enrich_gene", "ppi_source", "ppi_score", "term_id",
    "term_name", "group", "p_value", "p_adjust", "path_type"
  )
  paths <- list()

  if (nrow(ppi)) {
    directed_ppi <- rbind(
      data.frame(
        relation_gene = ppi$gene_a,
        enrich_gene = ppi$gene_b,
        ppi_source = ppi$ppi_source,
        ppi_score = ppi$ppi_score,
        stringsAsFactors = FALSE
      ),
      data.frame(
        relation_gene = ppi$gene_b,
        enrich_gene = ppi$gene_a,
        ppi_source = ppi$ppi_source,
        ppi_score = ppi$ppi_score,
        stringsAsFactors = FALSE
      )
    )
    ppi_paths <- merge(mpi, directed_ppi, by = "relation_gene", sort = FALSE)
    ppi_paths <- merge(
      ppi_paths,
      enrichment,
      by.x = "enrich_gene",
      by.y = "gene",
      sort = FALSE
    )
    if (nrow(ppi_paths)) {
      ppi_paths$path_type <- "ppi"
      paths[[length(paths) + 1L]] <- ppi_paths[common_columns]
    }
  }

  if (include_overlap) {
    overlap <- merge(
      mpi,
      enrichment,
      by.x = "relation_gene",
      by.y = "gene",
      sort = FALSE
    )
    if (nrow(overlap)) {
      overlap$enrich_gene <- overlap$relation_gene
      overlap$ppi_source <- "direct_membership"
      overlap$ppi_score <- 1
      overlap$path_type <- "gene_overlap"
      paths[[length(paths) + 1L]] <- overlap[common_columns]
    }
  }

  if (!length(paths)) return(.empty_path_table())
  out <- unique(do.call(rbind, paths))

  degree <- table(c(ppi$gene_a, ppi$gene_b))
  out$term_gene_degree <- as.integer(degree[out$enrich_gene])
  out$term_gene_degree[is.na(out$term_gene_degree)] <- 0L
  max_degree <- max(c(1L, as.integer(degree)))
  out$hub_factor <- switch(
    hub_mode,
    none = rep(1, nrow(out)),
    favor = 0.5 + 0.5 * log1p(out$term_gene_degree) / log1p(max_degree),
    penalize = 1 / sqrt(pmax(1, out$term_gene_degree))
  )
  out$hub_factor[out$path_type == "gene_overlap"] <- 1
  out$edge_weight <- out$mpi_weight * out$ppi_score * out$hub_factor
  out <- out[order(out$term_id, out$metabolite_name, -out$edge_weight), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.score_pairs <- function(metabolites, mpi, enrichment, paths, keep_zero) {
  term_table <- enrichment[!duplicated(enrichment$term_id),
    c("term_id", "term_name", "group", "p_value", "p_adjust"), drop = FALSE]
  relation_count <- aggregate(
    relation_gene ~ metabolite_key,
    mpi,
    function(x) length(unique(x))
  )
  names(relation_count)[names(relation_count) == "relation_gene"] <- "n_relation_genes"
  term_count <- aggregate(gene ~ term_id, enrichment, function(x) length(unique(x)))
  names(term_count)[names(term_count) == "gene"] <- "n_term_genes"
  term_table <- merge(term_table, term_count, by = "term_id", sort = FALSE)
  combinations <- .cartesian(
    metabolites[c("metabolite_key", "metabolite_id", "metabolite_name")],
    term_table
  )
  combinations <- merge(
    combinations,
    relation_count,
    by = "metabolite_key",
    sort = FALSE
  )

  rows <- lapply(seq_len(nrow(combinations)), function(i) {
    item <- combinations[i, , drop = FALSE]
    selected <- paths[
      paths$metabolite_key == item$metabolite_key & paths$term_id == item$term_id,
      , drop = FALSE
    ]
    edge_count <- nrow(selected)
    weighted_sum <- if (edge_count) sum(selected$edge_weight) else 0
    connected_relation <- if (edge_count) length(unique(selected$relation_gene)) else 0L
    connected_term <- if (edge_count) length(unique(selected$enrich_gene)) else 0L
    possible_pairs <- item$n_relation_genes * item$n_term_genes
    data.frame(
      item,
      edge_count = edge_count,
      weighted_sum = weighted_sum,
      normalized_score = if (possible_pairs > 0) weighted_sum / possible_pairs else NA_real_,
      mean_edge_weight = if (edge_count) weighted_sum / edge_count else 0,
      possible_pairs = possible_pairs,
      n_connected_relation_genes = connected_relation,
      n_connected_term_genes = connected_term,
      relation_coverage = connected_relation / item$n_relation_genes,
      term_coverage = connected_term / item$n_term_genes,
      stringsAsFactors = FALSE
    )
  })
  scores <- do.call(rbind, rows)
  if (!keep_zero) scores <- scores[scores$edge_count > 0, , drop = FALSE]
  scores$rank_within_term <- ave(
    scores$normalized_score,
    scores$term_id,
    FUN = .safe_rank
  )
  scores$global_rank <- .safe_rank(scores$normalized_score)
  scores <- scores[order(scores$term_id, scores$rank_within_term, scores$metabolite_name), , drop = FALSE]
  rownames(scores) <- NULL
  list(scores = scores, term_table = term_table)
}
