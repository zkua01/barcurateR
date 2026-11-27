# Compile QC flags and run the full QC pipeline

- `rb_compile_qc_flags()`: Combines logical check columns into a
  pipe-delimited `qc_flag` column.

- `rb_run_qc_pipeline()`: Runs the full QC pipeline in one call.

## Usage

``` r
rb_compile_qc_flags(data, exclude_flags = NULL, max_ambiguous_frac = NULL)

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

- exclude_flags:

  Character vector of flags that prevent `"pass"`. `NULL` (default) uses
  the hard-flag set defined in the function body.

- max_ambiguous_frac:

  Numeric scalar in `[0, 1]`, or `NULL` (default). Optional exclusion
  threshold on the proportion of ambiguous (IUPAC) bases among called
  bases (`frac_ambig`). When supplied, records with
  `frac_ambig > max_ambiguous_frac` receive the `ambiguous_excess` flag
  and are excluded from the passing set. When `NULL`, ambiguity never
  excludes a record; it is reported only through the informational
  `ambiguous_bases` note and the `n_ambiguous` / `frac_ambig` metadata
  columns. Requires the `frac_ambig` column to be present.

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

The input data frame `data` with two columns added:

- qc_flag:

  Character. Pipe-delimited excluding flags, or `"pass"` if none fired.
  This is the field consumed by
  [`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md).

- qc_note:

  Character. Pipe-delimited informational flags, or `NA_character_` if
  none fired.

## Details

Flags are partitioned into two classes. *Hard* flags (contaminant, numt,
internal_stop, frameshifted, sequence_short, divergent,
manual_review_needed, check_failed) prevent a record from being labelled
`"pass"`. *Soft* flags (sequence_gap, ambiguous_bases, frame_incomplete,
tree_error) describe properties or tooling limitations that do not
disqualify a record. Only excluding flags are written to `qc_flag`.

## See also

Other quality control:
[`rb_qc_checks`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md),
[`rb_qc_screening`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md)

## Examples

``` r
d <- data.frame(
  is_contaminant = FALSE, is_numt = FALSE,
  has_stop = FALSE, frameshifted = FALSE, has_short = FALSE,
  div_result = NA_character_, seq_type = "12S",
  n_gap_chars = 0L, n_ambiguous = 2L, frac_ambig = 0.01,
  frame_complete = NA, codon_check_run = NA
)
out <- rb_compile_qc_flags(d)
out$qc_flag   # "pass"
#> [1] "pass"
out$qc_note   # "ambiguous_bases"
#> [1] "ambiguous_bases"

d$is_contaminant <- TRUE
rb_compile_qc_flags(d)$qc_flag   # "contaminant"
#> [1] "contaminant"
```
