# Human and mouse MPI networks

`MPI_human_KEGG_Reactome_HumanGEM_BRENDA.csv` preserves the integrated human MPI evidence supplied in the local workbook. `MPI_mouse_Ensembl_orthologs.csv` is an Ensembl Compara orthology projection of that human evidence. Use `metabolite_kegg_id` and `gene_symbol` as the edge endpoints. `direction` is deliberately `not_inferred_from_source_workbook`: the workbook did not provide reaction-side information sufficient to label substrate/product direction. See `MPI_build_metadata.json` for counts, source hash, and limitations.
