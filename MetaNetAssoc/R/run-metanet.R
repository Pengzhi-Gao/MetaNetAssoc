#' Associate metabolites with biological terms through MPI and PPI networks
#'
#' `run_metanet()` evaluates direct PPI links between each metabolite's MPI genes
#' and each enriched term's genes. It retains every supporting path and reports
#' both raw counts and scores normalized by the theoretical number of gene pairs.
#'
#' @param metabolites Candidate metabolite IDs and/or names.
#' @param mpi Standardized MPI data from [prepare_mpi()].
#' @param ppi Standardized PPI data from [prepare_ppi()].
#' @param enrichment Standardized term-gene data from [prepare_enrichment()] or
#'   [run_enrichment()].
#' @param focus_genes Optional original gene set, used for QC reporting.
#' @param mpi_sources,ppi_sources Optional source names to retain.
#' @param include_overlap Treat an MPI gene that is itself a term member as a
#'   direct-membership path with PPI confidence one.
#' @param hub_mode `"none"` (default), `"favor"`, or `"penalize"`. Hub weighting
#'   is based on term-gene degree and always remains in the 0-1 range.
#' @param term_similarity_cutoff Optional Jaccard cutoff for greedy term reduction.
#' @param keep_zero Keep metabolite-term pairs with no supporting paths.
#' @return A [MetaNetResult-class] object.
#' @export
run_metanet <- function(
    metabolites = NULL,
    mpi = NULL,
    ppi = NULL,
    enrichment = NULL,
    focus_genes = NULL,
    mpi_sources = NULL,
    ppi_sources = NULL,
    include_overlap = TRUE,
    hub_mode = c("none", "favor", "penalize"),
    term_similarity_cutoff = NULL,
    keep_zero = TRUE,
    deg = NULL,
    species = c("human", "mouse"),
    mpi_database = NULL,
    string_data = NULL,
    string_score_threshold = 400) {
  call <- match.call()
  hub_mode <- match.arg(hub_mode)
  species <- match.arg(species)
  if (!is.null(deg)) {
    if (missing(metabolites) || is.null(metabolites)) {
      stop("Supply HMDB metabolite IDs in `metabolites` when using `deg`.", call. = FALSE)
    }
    if (is.null(ppi)) ppi <- build_ppi_network(deg, species, score_threshold = string_score_threshold,
      string_data = string_data)
    if (is.null(mpi)) mpi <- build_mpi_network(deg, metabolites, species, mpi_database = mpi_database)
    if (is.null(enrichment)) {
      enrichment <- run_deg_enrichment(deg, species)$enrichment
    }
    focus_genes <- .extract_deg(deg, "gene")
  }
  if (is.list(enrichment) && "enrichment" %in% names(enrichment)) enrichment <- enrichment$enrichment
  metabolites <- .clean_text(metabolites)
  metabolites <- unique(metabolites[!is.na(metabolites)])
  if (!length(metabolites)) stop("Supply at least one metabolite.", call. = FALSE)

  .assert_columns(
    mpi,
    c("metabolite_id", "metabolite_name", "relation_gene", "mpi_source", "mpi_weight"),
    "mpi"
  )
  .assert_columns(ppi, c("gene_a", "gene_b", "ppi_source", "ppi_score"), "ppi")
  .assert_columns(
    enrichment,
    c("term_id", "term_name", "gene", "group", "p_value", "p_adjust"),
    "enrichment"
  )

  mpi <- prepare_mpi(
    mpi,
    source = "mpi_source",
    weight = "mpi_weight",
    weight_scale = "unit"
  )
  if (!is.null(mpi_sources)) {
    mpi <- mpi[mpi$mpi_source %in% mpi_sources, , drop = FALSE]
  }
  if (!nrow(mpi)) stop("No MPI rows remained after source filtering.", call. = FALSE)

  if (!is.null(ppi_sources)) {
    ppi <- ppi[ppi$ppi_source %in% ppi_sources, , drop = FALSE]
  }
  ppi <- prepare_ppi(
    ppi,
    source = "ppi_source",
    score = "ppi_score",
    score_scale = "unit"
  )
  enrichment <- prepare_enrichment(
    enrichment,
    group = "group",
    p_value = "p_value",
    p_adjust = "p_adjust",
    fdr_cutoff = NULL,
    pvalue_cutoff = NULL
  )

  wanted_upper <- toupper(metabolites)
  wanted_lower <- tolower(metabolites)
  mpi_id <- toupper(mpi$metabolite_id)
  mpi_name <- tolower(mpi$metabolite_name)
  selected_rows <- (!is.na(mpi_id) & mpi_id %in% wanted_upper) |
    (!is.na(mpi_name) & mpi_name %in% wanted_lower)
  matched_input <- vapply(metabolites, function(value) {
    toupper(value) %in% mpi_id || tolower(value) %in% mpi_name
  }, logical(1))
  unmatched_metabolites <- metabolites[!matched_input]
  mpi <- mpi[selected_rows, , drop = FALSE]
  if (!nrow(mpi)) {
    stop("None of the requested metabolites was found in the MPI table.", call. = FALSE)
  }
  mpi <- .aggregate_mpi(mpi)
  selected_metabolites <- unique(mpi[c("metabolite_key", "metabolite_id", "metabolite_name")])

  if (!is.null(term_similarity_cutoff)) {
    reduction <- reduce_terms(enrichment, term_similarity_cutoff)
    enrichment <- reduction$enrichment
    term_reduction <- reduction$mapping
  } else {
    term_reduction <- data.frame(
      term_id = unique(enrichment$term_id),
      representative_term_id = unique(enrichment$term_id),
      jaccard = 1,
      retained = TRUE,
      stringsAsFactors = FALSE
    )
  }

  paths <- .build_paths(mpi, ppi, enrichment, include_overlap, hub_mode)
  scored <- .score_pairs(selected_metabolites, mpi, enrichment, paths, keep_zero)
  ppi_nodes <- unique(c(ppi$gene_a, ppi$gene_b))
  focus_genes <- unique(.clean_gene(focus_genes))
  focus_genes <- focus_genes[!is.na(focus_genes)]
  unmatched <- list(
    metabolites = unmatched_metabolites,
    mpi_genes_not_in_ppi = setdiff(unique(mpi$relation_gene), ppi_nodes),
    term_genes_not_in_ppi = setdiff(unique(enrichment$gene), ppi_nodes),
    focus_genes_without_retained_term = setdiff(focus_genes, unique(enrichment$gene))
  )
  qc <- data.frame(
    metric = c(
      "requested_metabolites", "matched_metabolites", "unmatched_metabolites",
      "mpi_relation_genes", "ppi_nodes", "retained_terms", "supported_paths",
      "zero_score_pairs"
    ),
    value = c(
      length(metabolites), nrow(selected_metabolites), length(unmatched_metabolites),
      length(unique(mpi$relation_gene)), length(ppi_nodes), nrow(scored$term_table),
      nrow(paths), sum(scored$scores$edge_count == 0)
    ),
    stringsAsFactors = FALSE
  )

  new(
    "MetaNetResult",
    call = call,
    parameters = list(
      mpi_sources = mpi_sources,
      ppi_sources = ppi_sources,
      include_overlap = include_overlap,
      hub_mode = hub_mode,
      term_similarity_cutoff = term_similarity_cutoff,
      keep_zero = keep_zero
    ),
    scores = scored$scores,
    edges = paths,
    term_table = scored$term_table,
    qc = qc,
    unmatched = unmatched,
    term_reduction = term_reduction
  )
}
