#' Run GO enrichment with clusterProfiler
#'
#' This is an optional convenience wrapper. The core network algorithm does not
#' require clusterProfiler and can use imported MetaScape or other enrichment
#' results prepared with [prepare_enrichment()].
#'
#' @param genes Gene symbols.
#' @param organism Human or mouse.
#' @param ontology GO ontology passed to `clusterProfiler::enrichGO`.
#' @param background Optional background gene symbols.
#' @param p_adjust_cutoff Benjamini-Hochberg adjusted p-value cutoff.
#' @param min_gene_set_size,max_gene_set_size Term size limits.
#' @return A standardized long enrichment table.
#' @export
run_enrichment <- function(
    genes,
    organism = c("human", "mouse"),
    ontology = "BP",
    background = NULL,
    p_adjust_cutoff = 0.05,
    min_gene_set_size = 10,
    max_gene_set_size = 500) {
  organism <- match.arg(organism)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install `clusterProfiler` to run enrichment, or import existing results.", call. = FALSE)
  }
  annotation_package <- if (organism == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
  if (!requireNamespace(annotation_package, quietly = TRUE)) {
    stop("Install `", annotation_package, "` to run enrichment for ", organism, ".", call. = FALSE)
  }
  orgdb <- getExportedValue(annotation_package, annotation_package)
  genes <- unique(.clean_gene(genes))
  genes <- genes[!is.na(genes)]
  background <- if (is.null(background)) NULL else unique(.clean_gene(background))
  if (!is.null(background)) background <- background[!is.na(background)]

  result <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = orgdb,
    keyType = "SYMBOL",
    ont = ontology,
    universe = background,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    readable = FALSE
  )
  result <- as.data.frame(result)
  if (!nrow(result)) {
    stop("No enriched terms were returned.", call. = FALSE)
  }
  result$ontology <- ontology
  out <- prepare_enrichment(
    result,
    term_id = "ID",
    term_name = "Description",
    gene = NULL,
    gene_list = "geneID",
    group = "ontology",
    p_value = "pvalue",
    p_adjust = "p.adjust",
    gene_sep = "/",
    fdr_cutoff = p_adjust_cutoff,
    pvalue_cutoff = NULL
  )
  attr(out, "raw_result") <- result
  out
}

#' Reduce redundant terms using member-gene Jaccard similarity
#'
#' Terms are considered from strongest to weakest adjusted p-value. A term whose
#' gene set is too similar to an already retained term is mapped to that retained
#' representative.
#'
#' @param enrichment A standardized long enrichment table.
#' @param jaccard_cutoff Similarity threshold in the interval 0-1.
#' @param within_group Only compare terms within the same group.
#' @return A list containing reduced `enrichment` and a `mapping` table.
#' @export
reduce_terms <- function(enrichment, jaccard_cutoff = 0.7, within_group = TRUE) {
  .assert_data_frame(enrichment, "enrichment")
  .assert_columns(
    enrichment,
    c("term_id", "term_name", "gene", "group", "p_value", "p_adjust"),
    "enrichment"
  )
  if (length(jaccard_cutoff) != 1L || !is.finite(jaccard_cutoff) ||
      jaccard_cutoff < 0 || jaccard_cutoff > 1) {
    stop("`jaccard_cutoff` must be between 0 and 1.", call. = FALSE)
  }

  term_info <- enrichment[!duplicated(enrichment$term_id),
    c("term_id", "term_name", "group", "p_value", "p_adjust"), drop = FALSE]
  priority <- ifelse(is.na(term_info$p_adjust), term_info$p_value, term_info$p_adjust)
  priority[is.na(priority)] <- Inf
  term_info <- term_info[order(priority, term_info$term_id), , drop = FALSE]
  gene_sets <- split(enrichment$gene, enrichment$term_id)
  kept <- character()
  mapping <- vector("list", nrow(term_info))

  for (i in seq_len(nrow(term_info))) {
    candidate <- term_info$term_id[i]
    comparable <- kept
    if (within_group && length(comparable)) {
      group_lookup <- setNames(term_info$group, term_info$term_id)
      comparable <- comparable[group_lookup[comparable] == term_info$group[i]]
    }
    representative <- candidate
    similarity <- 1
    if (length(comparable)) {
      similarities <- vapply(comparable, function(k) {
        a <- unique(gene_sets[[candidate]])
        b <- unique(gene_sets[[k]])
        length(intersect(a, b)) / length(union(a, b))
      }, numeric(1))
      best <- which.max(similarities)
      if (similarities[best] >= jaccard_cutoff) {
        representative <- comparable[best]
        similarity <- similarities[best]
      }
    }
    if (representative == candidate) kept <- c(kept, candidate)
    mapping[[i]] <- data.frame(
      term_id = candidate,
      representative_term_id = representative,
      jaccard = similarity,
      retained = representative == candidate,
      stringsAsFactors = FALSE
    )
  }

  mapping <- do.call(rbind, mapping)
  reduced <- enrichment[enrichment$term_id %in% kept, , drop = FALSE]
  rownames(reduced) <- NULL
  rownames(mapping) <- NULL
  list(enrichment = reduced, mapping = mapping)
}
