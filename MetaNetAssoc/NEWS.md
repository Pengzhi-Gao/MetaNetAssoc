# MetaNetAssoc 0.4.0

- Corrected STRING v12 API addressing (`version-12-0.string-db.org`) and added
  a system-curl fallback for Windows TLS failures; STRING API 0--1 and
  downloaded-file 0--1000 scores are now detected automatically.
- Added MPI matching diagnostics, available as
  `attr(mpi, "matching_diagnostics")`, and corrected the DEG/HMDB tutorial so
  its example inputs have valid MPI overlap.

## 0.3.0

- Added packaged human MPI evidence (KEGG, Reactome, Human-GEM, and BRENDA) and
  the Ensembl-Compara-projected mouse MPI network.
- Added `load_mpi_reference()` for loading either network with provenance
  metadata. Both networks use KEGG compound IDs; the mouse network remains an
  orthology projection rather than independent mouse curation.

## 0.2.0

- Added a DEG/HMDB-first workflow: `build_ppi_network()` retrieves or imports
  STRING PPI edges, and `build_mpi_network()` filters a traceable HMDB--gene
  evidence database by differential genes and metabolites.
- Added `run_deg_enrichment()` and `export_enrichment()` for clusterProfiler
  analysis and raw enrichment-table export.
- Added human/mouse selection, offline STRING input, a synthetic HMDB MPI
  reference example, and automatic construction of missing inputs in
  `run_metanet(deg = ...)`.

## 0.1.0

- Initial standard R package structure.
- Standardized MPI, PPI, and enrichment import helpers.
- Auditable direct PPI and gene-overlap paths.
- Raw a2, edge-weighted, normalized, and coverage scores.
- Optional hub favoring or penalization.
- Redundant-term reduction by member-gene Jaccard similarity.
- Ranking, heatmap, and four-layer network visualizations.
- Deterministic synthetic germ-cell example data and test suite.
