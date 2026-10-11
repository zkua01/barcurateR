# Taxonomy standardization, joining, and string formatting

Functions to standardize species names, join taxonomic lineage onto
reference data, and format taxonomy strings for export.

- `rb_build_taxonomy_string()`: Formats the seven taxonomic ranks into a
  semicolon-delimited string, optionally with rank prefixes.

- `rb_standardize_species()`: Maps alternative species names to valid
  names using a reference table, with an optional external taxonomy
  lookup callback as fallback.

- `rb_join_taxonomy()`: Joins a taxonomy table onto reference data,
  filling missing genus/family/order with a configurable label and
  applying constant higher-rank defaults.

## Usage

``` r
rb_build_taxonomy_string(data, style = c("plain", "rank_prefix"))

rb_standardize_species(
  names,
  reference_table,
  name_col = "alternative_clean",
  valid_col = "valid_name",
  taxonomy_lookup = NULL
)

rb_join_taxonomy(
  data,
  taxonomy_table,
  by = "species",
  defaults = list(),
  unassigned_label = "UNASSIGNED"
)
```

## Arguments

- data:

  A data frame containing at least the seven taxonomic ranks (`kingdom`
  through `species`) for `rb_build_taxonomy_string()`, or at least the
  column specified by `by` for `rb_join_taxonomy()`.

- style:

  Formatting style for `rb_build_taxonomy_string()`. One of `"plain"`
  (default) or `"rank_prefix"`.

- names:

  Character vector of species names to standardize.

- reference_table:

  A data frame mapping alternative names to valid names for
  `rb_standardize_species()`.

- name_col:

  Column in `reference_table` holding alternative/known names. Default
  `"alternative_clean"`.

- valid_col:

  Column in `reference_table` holding the valid name. Default
  `"valid_name"`.

- taxonomy_lookup:

  Optional callback `function(name) -> valid_name` used when
  `reference_table` has no match (e.g., an external catalog API). `NULL`
  by default.

- taxonomy_table:

  A data frame containing taxonomic lineage to join onto `data` for
  `rb_join_taxonomy()`. Must contain the column specified by `by`, and
  optionally columns for `genus`, `family`, `order`, and higher ranks.

- by:

  Column name to join on. Default `"species"`.

- defaults:

  Named list of constant higher-rank values applied to every row, e.g.,
  `list(kingdom = "Metazoa", phylum = "Chordata")`.

- unassigned_label:

  Label used when genus/family/order cannot be resolved. Default
  `"UNASSIGNED"`.

## Value

- `rb_build_taxonomy_string()` returns a character vector of formatted
  taxonomy strings, one per row.

- `rb_standardize_species()` returns a character vector of standardized
  species names, with `NA` for unmatched names.

- `rb_join_taxonomy()` returns the input data frame with taxonomic
  lineage columns joined.

## Details

`rb_build_taxonomy_string()` supports two styles:

|                 |                                                 |
|-----------------|-------------------------------------------------|
| Style           | Example output                                  |
| `"plain"`       | `Animalia;Chordata;Actinopterygii;...`          |
| `"rank_prefix"` | `k__Animalia;p__Chordata;c__Actinopterygii;...` |

Missing or empty values are replaced with `"unassigned"`.

`rb_join_taxonomy()` guards against silent Cartesian joins by checking
for duplicate keys in `taxonomy_table` before merging. If duplicates are
found, it stops with an informative error.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build taxonomy strings for DADA2 export
tax_strings <- rb_build_taxonomy_string(refs, style = "plain")

# Standardize species names
valid_names <- rb_standardize_species(
  raw_names,
  reference_table = name_mapping
)

# Join taxonomy onto curated data
final_data <- rb_join_taxonomy(
  qc_data,
  taxonomy_table = tax_table,
  defaults = list(kingdom = "Animalia", phylum = "Chordata")
)
} # }
```
