# Synthetic example data

These files form a deterministic, biologically plausible example of a germ-cell
multi-omics analysis. Metabolite and term IDs beginning with `SIM` or `GO:SIM`
are deliberately synthetic. Gene symbols are recognizable human-style symbols,
but all MPI/PPI edges and confidence scores are tutorial data rather than curated
evidence. Do not use them to support biological conclusions.

- `example_metabolites.csv`: candidate metabolites and the expected dominant
  signal used in tests.
- `example_focus_genes.csv`: a stage-focused gene set that could have produced
  the enrichment result.
- `example_mpi.csv`: metabolite-gene evidence from two mock sources.
- `example_ppi.csv`: weighted undirected interactions, including reciprocal rows
  to exercise deduplication.
- `example_enrichment.csv`: significant terms in long term-gene format.
# Packaged MPI references (version 0.3.0)

`mpi_human_kegg_reactome_humangem_brenda.csv` contains the integrated human
MPI evidence. `mpi_mouse_ensembl_orthologs.csv` is its Ensembl Compara
orthology projection to mouse. Both use KEGG compound IDs. See
`mpi_build_metadata.json` and `mpi_networks_README.md` for provenance and
limitations.

