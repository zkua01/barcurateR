# Compile QC flags and run the full QC pipeline

- `rb_compile_qc_flags()`: Combines logical check columns into a
  pipe-delimited `qc_flag` column.

- `rb_run_qc_pipeline()`: Runs the full QC pipeline in one call.

## Usage

``` r
rb_compile_qc_flags(
  data,
  checks = list(contaminant = "is_contaminant", numt = "is_numt", internal_stop =
    "has_stop", frameshifted = "frameshifted", sequence_short = "has_short", sequence_gap
    = "has_gaps"),
  extra_flag_fn = NULL
)

rb_run_qc_pipeline(
  data,
  blast_db,
  numt_fasta = NULL,
  species_col = "species",
  gene_col = "seq_type",
  sequence_col = "sequence",
  id_col = "unique_code",
  min_divergence_seqs = 5,
  parallel_workers = 1,
  coding_genes = c("COI"),
  rrna_genes = c("12S", "16S")
)
```

## Arguments

- data:

  A data frame containing sequences and check results.

- checks:

  Named list mapping flag names to logical column names.

- extra_flag_fn:

  Optional function returning additional flags.

- blast_db:

  Path to a BLAST database.

- numt_fasta:

  Path to a NUMT FASTA file.

- species_col, gene_col, sequence_col, id_col:

  Column names.

- min_divergence_seqs:

  Minimum sequences for divergence analysis.

- parallel_workers:

  Number of parallel workers.

- coding_genes, rrna_genes:

  Gene types for checks.

## Value

The input data frame with a `qc_flag` column added.

## See also

Other quality control:
[`rb_qc_checks`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md),
[`rb_qc_screening`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md)
