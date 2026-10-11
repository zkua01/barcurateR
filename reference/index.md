# Package index

## Database Curation & Standardization

Functions for standardizing, merging, and quality-controlling reference
sequences.

- [`rb_read_fasta()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_score_assignments()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_exact_hits()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_parse_blast_tabular()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_assign_blastn()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_assign_edna()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_build_species_matrix()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  : Assign eDNA sequences to reference taxa
- [`rb_default_marker_patterns()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md)
  [`rb_extract_markers()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md)
  [`rb_parse_source_table()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md)
  : Parse and standardize source sequence tables
- [`rb_combine_sources()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md)
  [`rb_resolve_ambiguous()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md)
  : Combine parsed sources and resolve ambiguous sequences
- [`rb_curate_reference()`](https://zkua01.github.io/barcurateR/reference/rb_curate_reference.md)
  : Curate a DNA reference sequence database
- [`rb_reverse_complement()`](https://zkua01.github.io/barcurateR/reference/rb_sequence_prep.md)
  [`rb_remove_primers()`](https://zkua01.github.io/barcurateR/reference/rb_sequence_prep.md)
  : Prepare and trim raw DNA sequences
- [`rb_detect_ambiguity()`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md)
  [`rb_write_ambiguity_csv()`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md)
  [`rb_read_ambiguity_csv()`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md)
  : Detect and archive ambiguous sequences

## Database I/O & Connection

Connect to, read from, and write to the SQLite database.

- [`rb_available_datasets()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md)
  [`rb_connect()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md)
  [`rb_disconnect()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md)
  [`rb_tables()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md)
  [`rb_read_table()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md)
  : Connect to and manage barcurateR SQLite databases
- [`rb_write_reference_table()`](https://zkua01.github.io/barcurateR/reference/rb_write_db.md)
  [`rb_write_qc_table()`](https://zkua01.github.io/barcurateR/reference/rb_write_db.md)
  [`rb_write_barcode_gap_table()`](https://zkua01.github.io/barcurateR/reference/rb_write_db.md)
  [`rb_write_ambiguity_table()`](https://zkua01.github.io/barcurateR/reference/rb_write_db.md)
  [`rb_build_reference_db()`](https://zkua01.github.io/barcurateR/reference/rb_write_db.md)
  : Write curated reference data to a SQLite database
- [`rb_yzfishdb_release()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  [`rb_default_yzfishdb_url()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  [`rb_db_url()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  [`rb_db_dir()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  [`rb_db_path()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  [`rb_db_available()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  [`rb_install_db()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)
  : Install and manage barcurateR reference databases
- [`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  [`rb_list_taxa()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  [`rb_filter_reference()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  [`rb_taxonomy()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  : Query and filter a curated reference database
- [`rb_export_fasta()`](https://zkua01.github.io/barcurateR/reference/rb_export.md)
  [`rb_export_blastn()`](https://zkua01.github.io/barcurateR/reference/rb_export.md)
  [`rb_export_dada2()`](https://zkua01.github.io/barcurateR/reference/rb_export.md)
  [`rb_export_qiime2()`](https://zkua01.github.io/barcurateR/reference/rb_export.md)
  [`rb_export_usearch()`](https://zkua01.github.io/barcurateR/reference/rb_export.md)
  : Export curated reference sequences to common metabarcoding formats

## Auditing & Diagnostics

Tools to inspect database health, marker coverage, and barcode gaps.

- [`rb_marker_coverage()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  [`rb_source_coverage()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  [`rb_qc_summary()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  [`rb_barcode_gap()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  [`rb_ambiguity()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  [`rb_read_barcode_gap()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  [`rb_read_ambiguity()`](https://zkua01.github.io/barcurateR/reference/rb_diagnostics.md)
  : Audit and inspect a curated reference database
- [`rb_check_codons()`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md)
  [`rb_check_rrna_integrity()`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md)
  [`rb_check_divergence()`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md)
  : Sequence integrity checks for quality control
- [`rb_compile_qc_flags()`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md)
  [`rb_run_qc_pipeline()`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md)
  : Compile QC flags and run the full QC pipeline
- [`rb_screen_contaminants()`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md)
  [`rb_screen_numts()`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md)
  : Screen reference sequences for contaminants and NUMTs
- [`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  [`rb_list_taxa()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  [`rb_filter_reference()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  [`rb_taxonomy()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
  : Query and filter a curated reference database
- [`rb_diagnostic_sites()`](https://zkua01.github.io/barcurateR/reference/rb_barcode_gap_analysis.md)
  [`rb_barcode_gap_species()`](https://zkua01.github.io/barcurateR/reference/rb_barcode_gap_analysis.md)
  [`rb_run_barcode_gap()`](https://zkua01.github.io/barcurateR/reference/rb_barcode_gap_analysis.md)
  : Barcode gap analysis for curated reference databases
- [`rb_build_gene_tree()`](https://zkua01.github.io/barcurateR/reference/rb_build_gene_tree.md)
  : Build a neighbor-joining gene tree for a marker

## Taxonomic Assignment

Assign eDNA ASVs to the curated database.

- [`rb_read_fasta()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_score_assignments()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_exact_hits()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_parse_blast_tabular()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_assign_blastn()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_assign_edna()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  [`rb_build_species_matrix()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
  : Assign eDNA sequences to reference taxa

## Utilities & Classification

Helper functions for taxonomy, ML classification, and contaminants.

- [`rb_build_contaminant_db()`](https://zkua01.github.io/barcurateR/reference/rb_build_contaminant_db.md)
  [`rb_build_blast_db()`](https://zkua01.github.io/barcurateR/reference/rb_build_contaminant_db.md)
  : Build contaminant reference databases for QC screening
- [`rb_sequence_features()`](https://zkua01.github.io/barcurateR/reference/rb_classification.md)
  [`rb_train_classifier()`](https://zkua01.github.io/barcurateR/reference/rb_classification.md)
  [`rb_classify_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_classification.md)
  : Machine learning classification of DNA sequences
- [`rb_build_taxonomy_string()`](https://zkua01.github.io/barcurateR/reference/rb_taxonomy_utils.md)
  [`rb_standardize_species()`](https://zkua01.github.io/barcurateR/reference/rb_taxonomy_utils.md)
  [`rb_join_taxonomy()`](https://zkua01.github.io/barcurateR/reference/rb_taxonomy_utils.md)
  : Taxonomy standardization, joining, and string formatting
