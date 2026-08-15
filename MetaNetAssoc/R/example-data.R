#' Load the biologically plausible synthetic example data
#'
#' The dataset mimics a germ-cell metabolomics/transcriptomics analysis. It is
#' synthetic and intended only for software testing and tutorials; associations
#' must not be interpreted as curated biological evidence.
#'
#' @param path Optional path to an `extdata` directory. Normally omitted.
#' @return A named list containing candidate metabolites, focus genes, MPI, PPI,
#'   and enrichment tables.
#' @export
load_example_data <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", package = "MetaNetAssoc")
  }
  if (!nzchar(path) || !dir.exists(path)) {
    stop("Example data directory was not found. Install the package first or supply `path`.", call. = FALSE)
  }
  read_example <- function(file) {
    utils::read.csv(
      file.path(path, file),
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA")
    )
  }
  metabolites <- read_example("example_metabolites.csv")
  focus_genes <- read_example("example_focus_genes.csv")$gene
  mpi_raw <- read_example("example_mpi.csv")
  ppi_raw <- read_example("example_ppi.csv")
  enrichment_raw <- read_example("example_enrichment.csv")
  hmdb_mpi_reference <- read_example("example_hmdb_mpi_reference.csv")

  list(
    metabolites = metabolites,
    hmdb_metabolites = unique(hmdb_mpi_reference[c("hmdb_id", "metabolite_name")]),
    hmdb_mpi_reference = hmdb_mpi_reference,
    focus_genes = focus_genes,
    mpi = prepare_mpi(
      mpi_raw,
      source = "mpi_source",
      weight = "mpi_weight",
      weight_scale = "unit"
    ),
    ppi = ppi_raw,
    enrichment = prepare_enrichment(
      enrichment_raw,
      group = "group",
      p_value = "p_value",
      p_adjust = "p_adjust",
      fdr_cutoff = 0.05,
      pvalue_cutoff = NULL
    )
  )
}

#' Load a packaged human or mouse MPI reference network
#'
#' The human network integrates the project-provided KEGG, Reactome, Human-GEM,
#' and BRENDA evidence. The mouse network projects those human edges through
#' Ensembl Compara orthology. Metabolites are represented by KEGG compound IDs,
#' not HMDB IDs; use a documented HMDB--KEGG cross-reference when an HMDB-first
#' workflow is required.
#'
#' @param species Either `"human"` or `"mouse"`.
#' @param include_metadata Return a list containing both the network and its
#'   build metadata. Otherwise return only the network table.
#' @return A data frame, or a list with `network` and `metadata`.
#' @export
load_mpi_reference <- function(species = c("human", "mouse"), include_metadata = FALSE) {
  species <- match.arg(species)
  file <- if (species == "human") {
    "mpi_human_kegg_reactome_humangem_brenda.csv"
  } else {
    "mpi_mouse_ensembl_orthologs.csv"
  }
  path <- system.file("extdata", file, package = "MetaNetAssoc")
  if (!nzchar(path)) stop("Packaged MPI reference was not found. Reinstall MetaNetAssoc.", call. = FALSE)
  network <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!include_metadata) return(network)
  meta_path <- system.file("extdata", "mpi_build_metadata.json", package = "MetaNetAssoc")
  metadata <- if (requireNamespace("jsonlite", quietly = TRUE) && nzchar(meta_path)) {
    jsonlite::fromJSON(meta_path, simplifyVector = TRUE)
  } else {
    list(metadata_file = meta_path, note = "Install jsonlite to parse metadata automatically.")
  }
  list(network = network, metadata = metadata)
}
