# Combine parsed sources and resolve ambiguous sequences

Functions to combine multiple standardized source tables and resolve
sequences that are assigned to multiple species.

- `rb_combine_sources()`: Combines parsed source tables and removes
  duplicated species-sequence combinations.

- `rb_resolve_ambiguous()`: Detects and resolves sequences assigned to
  multiple species.

## Usage

``` r
rb_combine_sources(
  source_list,
  species_col = "species",
  sequence_col = "sequence"
)

rb_resolve_ambiguous(
  data,
  species_col = "species",
  sequence_col = "sequence",
  resolution_table = NULL,
  on_unresolved = "stop",
  flag_output_path = NULL
)
```

## Arguments

- source_list:

  List of standardized source data frames.

- species_col:

  Name of the species column. Default `"species"`.

- sequence_col:

  Name of the sequence column. Default `"sequence"`.

- data:

  Standardized data frame to check for ambiguous sequences.

- resolution_table:

  Optional data frame resolving ambiguous sequences. Must contain an
  `action` column with values `"keep"`, `"combine"`, or `"drop"`.

- on_unresolved:

  Action to take when ambiguous sequences are found and no
  `resolution_table` is supplied. One of `"stop"`, `"warn"`, or
  `"drop"`. Default `"stop"`.

- flag_output_path:

  Optional path to write ambiguous rows for manual review before
  resolving.

## Value

- `rb_combine_sources()` returns a combined data frame with duplicated
  species-sequence combinations removed.

- `rb_resolve_ambiguous()` returns the input data frame with ambiguous
  rows removed or replaced according to the `resolution_table` (or
  `on_unresolved` rules if no table is provided).

## Details

Ambiguous sequences are sequences that appear under more than one
species. These must be resolved before QC because downstream QC and
curation assume that each sequence is assigned to a single species.

If a `resolution_table` is supplied, it should contain one row per
ambiguous sequence and an `action` column:

- `"keep"`: keep the specified resolution.

- `"combine"`: combine records across species, pipe-joining values that
  differ across the ambiguous group.

- `"drop"`: remove the ambiguous sequence.

If no `resolution_table` is supplied, `on_unresolved` controls the
behavior:

- `"stop"`: stop with an error.

- `"warn"`: warn and return the data unchanged.

- `"drop"`: warn and remove ambiguous rows.

## See also

Other parsing:
[`rb_ambiguity_utils`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md),
[`rb_parse_sources`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md),
[`rb_sequence_prep`](https://zkua01.github.io/barcurateR/reference/rb_sequence_prep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
combined <- rb_combine_sources(
  list(ncbi_parsed, bold_parsed),
  species_col = "species",
  sequence_col = "sequence"
)

resolved <- rb_resolve_ambiguous(
  combined,
  on_unresolved = "drop"
)
} # }
```
