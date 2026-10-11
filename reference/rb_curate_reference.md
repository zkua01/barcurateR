# Curate a DNA reference sequence database

`rb_curate_reference()` is the main user-facing wrapper that
orchestrates the full reference database curation pipeline:

1.  First-pass quality control (contaminants, NUMTs, codons, rRNA).

2.  Optional machine-learning classification of unknown sequences.

3.  Second-pass quality control using updated marker labels.

4.  Optional taxonomic lineage joining and validation.

5.  Optional barcode gap analysis.

6.  Optional SQLite database export.

## Usage

``` r
rb_curate_reference(
  data,
  blast_db = NULL,
  numt_fasta = NULL,
  run_classifier = FALSE,
  run_divergence = FALSE,
  run_barcode_gap = FALSE,
  taxonomy_table = NULL,
  ambiguity_data = NULL,
  db_path = NULL,
  confidence_threshold = 0.8,
  min_length = 100,
  classifier_trees = 500,
  barcode_gap_markers = NULL,
  barcode_gap_min_seqs = 5,
  barcode_gap_parallel = FALSE,
  reference_table = "reference_final",
  qc_table = "qc_reference",
  barcode_gap_table = "barcode_gap_metrics",
  ambiguity_table = "ambiguous_sequences",
  ambiguity_path = "ambiguous_sequences.csv",
  archive_ambiguity = TRUE,
  check_ambiguity = TRUE,
  on_ambiguous = c("stop", "warn", "ignore"),
  require_taxonomy = TRUE,
  mafft = "mafft",
  fasttree = "FastTree"
)
```

## Arguments

- data:

  A standardized data frame containing at least `species`, `sequence`,
  and `seq_type`. A `source` column is also recommended.

- blast_db:

  Path prefix to a pre-built BLAST database used for contaminant
  screening.

- numt_fasta:

  Path to a FASTA file of known NUMT sequences.

- run_classifier:

  Logical. If `TRUE`, trains a Random Forest classifier on known marker
  types and predicts marker types for sequences labeled `"other"` or
  `"other_unknown"`.

- run_divergence:

  Logical. If `TRUE`, runs phylogenetic divergence checks. Requires
  `mafft` and `FastTree` to be installed.

- run_barcode_gap:

  Logical. If `TRUE`, runs barcode gap analysis after QC and taxonomy
  joining.

- taxonomy_table:

  Optional data frame used to join taxonomic lineage onto the curated
  sequences. Must contain `species` and the seven taxonomic ranks
  (kingdom through species).

- ambiguity_data:

  Optional data frame of ambiguity records to archive in the SQLite
  database.

- db_path:

  Optional path to write the final SQLite database.

- confidence_threshold:

  Minimum prediction probability required for the ML classifier to
  reclassify an unknown sequence (only used if `run_classifier = TRUE`).

- min_length:

  Minimum sequence length required for ML classification (only used if
  `run_classifier = TRUE`).

- classifier_trees:

  Number of trees used by the Random Forest classifier (only used if
  `run_classifier = TRUE`).

- barcode_gap_markers:

  Optional character vector of markers to include in barcode gap
  analysis. If `NULL`, all markers except `"other"` are used (only used
  if `run_barcode_gap = TRUE`).

- barcode_gap_min_seqs:

  Minimum number of sequences per species required for barcode gap
  analysis (only used if `run_barcode_gap = TRUE`).

- barcode_gap_parallel:

  Logical. If `TRUE`, runs barcode gap analysis in parallel where
  supported (only used if `run_barcode_gap = TRUE`).

- reference_table:

  Name of the final reference table written to the SQLite database.
  Default `"reference_final"`.

- qc_table:

  Name of the QC table. Default `"qc_reference"`.

- barcode_gap_table:

  Name of the barcode gap table. Default `"barcode_gap_metrics"` (only
  used if `run_barcode_gap = TRUE`).

- ambiguity_table:

  Name of the ambiguity table. Default `"ambiguous_sequences"`.

- ambiguity_path:

  Path to a CSV file containing ambiguity records. If
  `ambiguity_data = NULL`, `archive_ambiguity = TRUE`, and this file
  exists, it will be read and archived.

- archive_ambiguity:

  Logical. If `TRUE`, attempts to archive ambiguity records in the
  SQLite database.

- check_ambiguity:

  Logical. If `TRUE`, checks the input data for sequences assigned to
  multiple species before running QC.

- on_ambiguous:

  Action to take if ambiguous input sequences are found. One of
  `"stop"`, `"warn"`, or `"ignore"`.

- require_taxonomy:

  Logical. If `TRUE`, requires the final curated table to contain full
  kingdom-to-species taxonomy.

- mafft:

  Path to the `mafft` executable (only used if `run_divergence = TRUE`).

- fasttree:

  Path to the `FastTree` executable (only used if
  `run_divergence = TRUE`).

## Value

An invisible list containing:

- qc1:

  Data frame after first-pass QC.

- classifier:

  Trained ML classifier object, or `NULL`.

- qc2:

  Data frame after second-pass QC.

- final_data:

  Final curated data frame.

- barcode_gap:

  Barcode gap metrics, or `NULL`.

- db_path:

  Path to the SQLite database, or `NULL`.

## Details

If `check_ambiguity = TRUE`, the wrapper checks whether the same
sequence is assigned to multiple species before QC begins. Such
ambiguities should normally be resolved during standardization using
[`rb_resolve_ambiguous()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md).

## See also

[`rb_run_qc_pipeline()`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md),
[`rb_build_reference_db()`](https://zkua01.github.io/barcurateR/reference/rb_write_db.md)

## Examples
