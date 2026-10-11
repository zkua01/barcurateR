# Barcode gap analysis for curated reference databases

Functions to assess whether reference sequences for a species form a
distinguishable barcode gap relative to sequences from other species.

- `rb_run_barcode_gap()`: Runs barcode gap analysis across one or more
  markers. Aligns sequences with
  [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html),
  then computes per-species distance metrics.

- `rb_barcode_gap_species()`: Computes barcode gap metrics for a single
  species within an aligned marker dataset.

- `rb_diagnostic_sites()`: Counts fixed diagnostic nucleotide positions
  that distinguish focal sequences from non-focal sequences.

## Usage

``` r
rb_diagnostic_sites(focal_seqs, other_seqs, min_freq = 0.9, max_other = 500)

rb_barcode_gap_species(
  species_name,
  data,
  aligned_seqs,
  seq_len,
  species_col = "species",
  genus_col = "genus",
  marker = NA_character_,
  min_seqs = 2,
  max_other = 2000,
  max_congeners = 500
)

rb_run_barcode_gap(
  data,
  markers = c("COI", "12S", "16S", "genome"),
  species_col = "species",
  genus_col = NULL,
  gene_col = "seq_type",
  sequence_col = "sequence",
  qc_flag_col = "qc_flag",
  min_seqs = 5,
  parallel = TRUE,
  n_workers = NULL,
  remove_gaps = TRUE
)
```

## Arguments

- focal_seqs:

  Character vector of focal sequences for `rb_diagnostic_sites()`.

- other_seqs:

  Character vector of non-focal sequences for `rb_diagnostic_sites()`.

- min_freq:

  Minimum focal base frequency required for a position to be considered
  diagnostic. Default `0.9`.

- max_other:

  Maximum number of non-focal sequences to use when computing global
  inter-specific distances. Default `500` in `rb_diagnostic_sites()` and
  `2000` in `rb_barcode_gap_species()`.

- species_name:

  Species to analyze in `rb_barcode_gap_species()`.

- data:

  A standardized reference data frame.

- aligned_seqs:

  Named character vector of aligned sequences.

- seq_len:

  Aligned sequence length used to normalize distances.

- species_col:

  Name of the species column. Default `"species"`.

- genus_col:

  Name of the genus column. If `NULL`, genus is derived from the first
  token of the species name.

- marker:

  Marker name associated with the aligned sequences.

- min_seqs:

  Minimum number of sequences required per species.

- max_congeners:

  Maximum number of congeneric sequences to use when computing
  congeneric distances. Default `500`.

- markers:

  Character vector of marker names to analyze.

- gene_col:

  Name of the marker/gene column. Default `"seq_type"`.

- sequence_col:

  Name of the sequence column. Default `"sequence"`.

- qc_flag_col:

  Name of the QC flag column. Default `"qc_flag"`. If present, only
  sequences with `qc_flag == "pass"` are analyzed.

- parallel:

  Logical. If `TRUE`, uses a parallel cluster where possible.

- n_workers:

  Optional number of parallel workers. If `NULL`, defaults to
  `parallel::detectCores() - 2`.

- remove_gaps:

  Logical. If `TRUE`, removes `"-"` characters from sequences before
  analysis. Default `TRUE`.

## Value

- `rb_run_barcode_gap()` returns a data frame of barcode gap metrics,
  with one row per species per marker. Returns an empty data frame if no
  suitable data are available.

- `rb_barcode_gap_species()` returns a one-row data frame of barcode gap
  metrics, or `NULL` if the species does not meet the minimum sequence
  threshold or has no non-focal comparison sequences.

- `rb_diagnostic_sites()` returns an integer count of diagnostic sites.

## Details

`rb_run_barcode_gap()` expects a standardized reference data frame with
at least species, marker/gene, and sequence columns. If a `qc_flag`
column is present, only sequences with `qc_flag == "pass"` are used.

For each marker, the function:

1.  Filters to the requested marker.

2.  Optionally removes alignment gap characters (`"-"`).

3.  Filters species with at least `min_seqs` sequences.

4.  Aligns sequences using
    [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html).

5.  Computes barcode gap metrics for each species.

The returned data frame includes intra- and inter-specific distance
summaries, gap metrics, congeneric comparisons where available,
diagnostic site counts, sampling adequacy, taxonomic resolution,
assignment risk, and a recommended similarity threshold.

If `parallel = TRUE`, a parallel cluster is used for the per-species
calculations. If cluster setup or execution fails, the function
automatically falls back to serial execution for that marker.

Distance calculations use `rb_string_dist()`, which prefers
[`pwalign::stringDist()`](https://rdrr.io/pkg/pwalign/man/stringDist.html)
and falls back to `Biostrings::stringDist()` on older Biostrings
versions.

When `parallel = TRUE`, progress bars are displayed via the `pbapply`
package.

## Output columns

The data frame returned by `rb_run_barcode_gap()` and
`rb_barcode_gap_species()` includes:

- `species`, `marker`, `n_sequences`

- `max_intra`, `min_inter`

- `gap_exists`, `gap_width`

- `p_overlap`, `q_overlap`

- `median_intra`, `median_inter`

- `intra_95`, `inter_5`, `percentile_gap`, `robust_gap`

- `n_congeners_in_db`, `min_congeneric`, `congeneric_gap`

- `congeneric_p_overlap`, `congeneric_q_overlap`

- `n_diagnostic_sites`

- `sampling_adequacy`, `taxonomic_resolution`

- `gap_may_be_undersampling_artifact`

- `error_msg`

- `assignment_risk`, `recommended_threshold`

## Examples

``` r
if (FALSE) { # \dontrun{
gap_results <- rb_run_barcode_gap(
  final_data,
  markers = c("COI", "12S"),
  min_seqs = 5,
  parallel = FALSE
)

head(gap_results)
} # }
```
