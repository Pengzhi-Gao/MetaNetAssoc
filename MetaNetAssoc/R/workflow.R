#' Build a STRING PPI network from differential genes
#'
#' Queries STRING's public API and returns an undirected, deduplicated PPI table.
#' A supplied `string_data` table makes the function fully offline and is useful
#' for reproducible analyses or restricted networks.
#'
#' @param deg A character vector of gene symbols or a data frame containing a
#'   gene-symbol column.
#' @param species Either `"human"` or `"mouse"`.
#' @param gene_column Column containing gene symbols when `deg` is a data frame.
#' @param score_threshold Minimum STRING combined score, on the 0--1000 scale.
#' @param string_data Optional STRING-like edge table with preferred names
#'   `preferredName_A`, `preferredName_B`, and `score`.
#' @param string_version STRING API version, e.g. `"12-0"` or `"12.0"`.
#' @return A standardized PPI data frame with STRING scores converted to 0--1.
#' @export
build_ppi_network <- function(deg, species = c("human", "mouse"),
                              gene_column = "gene", score_threshold = 400,
                              string_data = NULL, string_version = "12-0") {
  species <- match.arg(species)
  genes <- .extract_deg(deg, gene_column)
  if (!length(genes)) stop("No valid differential genes were supplied.", call. = FALSE)
  if (!is.numeric(score_threshold) || length(score_threshold) != 1L ||
      score_threshold < 0 || score_threshold > 1000) {
    stop("`score_threshold` must be between 0 and 1000.", call. = FALSE)
  }
  if (is.null(string_data)) {
    taxon <- if (species == "human") 9606L else 10090L
    ids <- utils::URLencode(paste(genes, collapse = "\r"), reserved = TRUE)
    version <- gsub("[._]", "-", as.character(string_version))
    endpoint <- paste0("https://version-", version,
      ".string-db.org/api/tsv/network?identifiers=", ids,
      "&species=", taxon, "&required_score=", as.integer(score_threshold),
      "&caller_identity=MetaNetAssoc")
    string_data <- .read_string_network(endpoint)
  }
  .assert_data_frame(string_data, "string_data")
  a <- .first_existing(string_data, c("preferredName_A", "gene_a", "protein1"))
  b <- .first_existing(string_data, c("preferredName_B", "gene_b", "protein2"))
  score <- .first_existing(string_data, c("score", "combined_score", "ppi_score"))
  if (is.null(a) || is.null(b) || is.null(score)) {
    stop("STRING data needs endpoint columns and a score column.", call. = FALSE)
  }
  raw <- data.frame(gene_a = string_data[[a]], gene_b = string_data[[b]],
    ppi_source = "STRING", ppi_score = string_data[[score]], stringsAsFactors = FALSE)
  score_values <- suppressWarnings(as.numeric(raw$ppi_score))
  # STRING's TSV API currently returns combined scores on 0--1, while exported
  # STRING files may use 0--1000. The user-facing threshold remains 0--1000.
  observed_max <- suppressWarnings(max(score_values, na.rm = TRUE))
  cutoff <- if (is.finite(observed_max) && observed_max <= 1) score_threshold / 1000 else score_threshold
  raw <- raw[is.finite(score_values) & score_values >= cutoff, , drop = FALSE]
  prepare_ppi(raw, source = "ppi_source", score = "ppi_score", score_scale = "auto")
}

#' Build an MPI network for differential genes and HMDB metabolites
#'
#' HMDB identifiers are used as stable metabolite keys. Because HMDB does not
#' provide a general public metabolite--gene bulk API, the biological evidence is
#' supplied through `mpi_database`; this keeps evidence provenance explicit.
#'
#' @param deg Differential genes as a vector or data frame.
#' @param metabolites HMDB IDs, e.g. `"HMDB0000122"`, or a data frame containing
#'   an HMDB ID column.
#' @param species `"human"` or `"mouse"`.
#' @param mpi_database Evidence table containing HMDB ID and gene columns. If
#'   omitted, the included synthetic reference data are used for demonstration.
#' @param gene_column,hmdb_column Column names for data-frame inputs.
#' @return A standardized MPI table restricted to the requested genes and HMDB IDs.
#' @export
build_mpi_network <- function(deg, metabolites, species = c("human", "mouse"),
                              mpi_database = NULL, gene_column = "gene",
                              hmdb_column = "hmdb_id") {
  species <- match.arg(species)
  genes <- .extract_deg(deg, gene_column)
  ids <- .extract_hmdb(metabolites, hmdb_column)
  if (is.null(mpi_database)) {
    mpi_database <- load_example_data()$mpi
    warning("Using the bundled synthetic MPI reference. Supply `mpi_database` for biological analysis.", call. = FALSE)
  }
  .assert_data_frame(mpi_database, "mpi_database")
  id_col <- .first_existing(mpi_database, c("hmdb_id", "metabolite_id"))
  gene_col <- .first_existing(mpi_database, c("relation_gene", "gene", "gene_id"))
  if (is.null(id_col) || is.null(gene_col)) stop("`mpi_database` needs HMDB ID and gene columns.", call. = FALSE)
  database_ids <- toupper(.clean_text(mpi_database[[id_col]]))
  database_genes <- .clean_gene(mpi_database[[gene_col]])
  matched_ids <- intersect(ids, unique(database_ids[!is.na(database_ids)]))
  matched_genes <- intersect(genes, unique(database_genes[!is.na(database_genes)]))
  rows <- database_ids %in% ids & database_genes %in% genes
  out <- mpi_database[rows, , drop = FALSE]
  diagnostics <- data.frame(
    metric = c("input_differential_genes", "input_hmdb_ids", "deg_symbols_in_mpi_database",
      "hmdb_ids_in_mpi_database", "matched_mpi_edges"),
    value = c(length(genes), length(ids), length(matched_genes), length(matched_ids), nrow(out)),
    stringsAsFactors = FALSE
  )
  if (!nrow(out)) {
    gene_example <- paste(head(matched_genes, 5), collapse = ", ")
    id_example <- paste(head(matched_ids, 5), collapse = ", ")
    warning(
      "No MPI edges matched the supplied DEG and HMDB IDs. Diagnostics: ",
      length(matched_genes), "/", length(genes), " DEG symbols occur in `mpi_database`",
      if (length(matched_genes)) paste0(" (e.g. ", gene_example, ")") else "",
      "; ", length(matched_ids), "/", length(ids), " HMDB IDs occur in `mpi_database`",
      if (length(matched_ids)) paste0(" (e.g. ", id_example, ")") else "",
      "; joint MPI edges: 0. Check that the same metabolite is linked to at least one supplied DEG.",
      call. = FALSE
    )
  }
  out <- prepare_mpi(out, metabolite_id = id_col,
    metabolite_name = if ("metabolite_name" %in% names(out)) "metabolite_name" else NULL,
    gene = gene_col,
    source = .first_existing(out, c("mpi_source", "source")),
    weight = .first_existing(out, c("mpi_weight", "weight", "score")))
  attr(out, "matching_diagnostics") <- diagnostics
  out
}

#' Run DEG enrichment and optionally export clusterProfiler's result table
#'
#' @param deg Differential genes as a vector or data frame.
#' @param species `"human"` or `"mouse"`.
#' @param export_file Optional CSV path for the unmodified clusterProfiler table.
#' @param ... Passed to [run_enrichment()].
#' @return A list with `raw` clusterProfiler-style table and `enrichment` long table.
#' @export
run_deg_enrichment <- function(deg, species = c("human", "mouse"),
                               gene_column = "gene", export_file = NULL, ...) {
  species <- match.arg(species)
  genes <- .extract_deg(deg, gene_column)
  enrichment <- run_enrichment(genes, organism = species, ...)
  raw <- attr(enrichment, "raw_result")
  if (is.null(raw)) raw <- enrichment
  if (!is.null(export_file)) export_enrichment(raw, export_file)
  list(raw = raw, enrichment = enrichment)
}

#' Export an enrichment result table as CSV
#' @param enrichment A data frame, or the list returned by [run_deg_enrichment()].
#' @param file Output CSV path.
#' @return The normalized output path, invisibly.
#' @export
export_enrichment <- function(enrichment, file) {
  if (is.list(enrichment) && !is.data.frame(enrichment) && "raw" %in% names(enrichment)) enrichment <- enrichment$raw
  .assert_data_frame(enrichment, "enrichment")
  utils::write.csv(enrichment, file, row.names = FALSE, na = "")
  invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
}

.extract_deg <- function(deg, gene_column) {
  x <- if (is.data.frame(deg)) {
    .assert_columns(deg, gene_column, "deg"); deg[[gene_column]]
  } else deg
  x <- unique(.clean_gene(x)); x[!is.na(x)]
}
.extract_hmdb <- function(metabolites, hmdb_column) {
  x <- if (is.data.frame(metabolites)) {
    .assert_columns(metabolites, hmdb_column, "metabolites"); metabolites[[hmdb_column]]
  } else metabolites
  x <- unique(toupper(.clean_text(x))); x <- x[!is.na(x)]
  if (!length(x) || any(!grepl("^HMDB[0-9]+$", x))) stop("`metabolites` must contain HMDB IDs such as HMDB0000122.", call. = FALSE)
  x
}
.first_existing <- function(x, candidates) { hit <- intersect(candidates, names(x)); if (length(hit)) hit[1] else NULL }

.read_string_network <- function(endpoint) {
  first_error <- NULL
  out <- tryCatch(
    utils::read.delim(endpoint, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { first_error <<- e; NULL }
  )
  if (!is.null(out)) return(out)

  # On Windows, R's TLS stack can differ from the system curl executable.
  curl <- Sys.which("curl")
  if (nzchar(curl)) {
    file <- tempfile(fileext = ".tsv")
    on.exit(unlink(file), add = TRUE)
    status <- tryCatch(system2(curl, c("--fail", "--location", "--silent", "--show-error",
      "--retry", "2", "--output", shQuote(file), shQuote(endpoint))),
      error = function(e) 1L)
    if (identical(as.integer(status), 0L) && file.exists(file) && file.info(file)$size > 0) {
      return(utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE))
    }
  }
  stop(
    "STRING query failed at ", endpoint, ". R reported: ", conditionMessage(first_error),
    ". Check network/TLS settings, or download STRING edges and supply `string_data` for offline use.",
    call. = FALSE
  )
}
