# Detect and archive ambiguous sequences

Utility functions for detecting sequences assigned to multiple species
and archiving ambiguity records for downstream auditing.

- `rb_detect_ambiguity()`: Flags sequences that are assigned to more
  than one species.

- `rb_write_ambiguity_csv()`: Writes an ambiguity table to CSV.

- `rb_read_ambiguity_csv()`: Reads an ambiguity table from CSV.

## Usage

``` r
rb_detect_ambiguity(data, species_col = "species", sequence_col = "sequence")

rb_write_ambiguity_csv(data, path = "ambiguous_sequences.csv")

rb_read_ambiguity_csv(path = "ambiguous_sequences.csv")
```

## Arguments

- data:

  A standardized data frame containing at least the columns specified by
  `species_col` and `sequence_col`.

- species_col:

  Name of the species column. Default `"species"`.

- sequence_col:

  Name of the sequence column. Default `"sequence"`.

- path:

  Path to the ambiguity CSV file. Default `"ambiguous_sequences.csv"`.

## Value

- `rb_detect_ambiguity()` returns the input data frame with an
  `is_ambiguous` logical column added.

- `rb_write_ambiguity_csv()` returns the output path invisibly.

- `rb_read_ambiguity_csv()` returns a data frame.

## Details

Ambiguous sequences are sequences that occur under more than one species
name. These should normally be resolved before QC using
[`rb_resolve_ambiguous()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md)
or a manual resolution table.

The ambiguity table is an auditing artifact. It is not used directly by
the QC pipeline, but it can be archived in the final SQLite database by
[`rb_curate_reference()`](https://zkua01.github.io/barcurateR/reference/rb_curate_reference.md)
using the `ambiguity_data` or `ambiguity_path` arguments.

## See also

Other parsing:
[`rb_combine_resolve`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md),
[`rb_parse_sources`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md),
[`rb_sequence_prep`](https://zkua01.github.io/barcurateR/reference/rb_sequence_prep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
combined <- rb_combine_sources(list(ncbi_parsed, bold_parsed))

combined <- rb_detect_ambiguity(combined)

ambiguities <- combined[combined$is_ambiguous, ]

if (nrow(ambiguities) > 0) {
  rb_write_ambiguity_csv(ambiguities, "ambiguous_sequences.csv")
}
} # }
```
