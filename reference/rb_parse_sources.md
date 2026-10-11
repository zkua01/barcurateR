# Parse and standardize source sequence tables

Functions to extract marker information from sequence descriptions and
parse raw source tables into the standardized barcurateR schema.

- `rb_default_marker_patterns()`: Returns the default regular
  expressions used to detect common mitochondrial and ribosomal markers.

- `rb_extract_markers()`: Extracts marker type and completeness from a
  free-text sequence description.

- `rb_parse_source_table()`: Parses a raw source table into the
  standardized schema used by downstream QC and curation functions.

## Usage

``` r
rb_default_marker_patterns()

rb_extract_markers(
  description,
  marker_patterns = rb_default_marker_patterns(),
  genome_pattern = "complete genome|mitochondrion, complete genome|mitogenome",
  multi_gene_pattern = paste0("genes for|gene for|and|;|,.*and|partial.*complete|",
    "complete.*partial|complete.*complete|rrna.*rrna|", "rrna.*trna|trna.*rrna"),
  other_gene_pattern = "trna|transfer rna|nd[0-9]|atp[0-9]|nad[0-9]"
)

rb_parse_source_table(
  raw,
  source_name,
  column_map,
  marker_fn = rb_extract_markers,
  description_col = NULL,
  min_length = 1
)
```

## Arguments

- description:

  Free-text sequence description, for example a GenBank definition line.

- marker_patterns:

  Named list mapping marker names to detection regular expressions.

- genome_pattern:

  Regular expression identifying complete genome records.

- multi_gene_pattern:

  Regular expression identifying multi-gene records.

- other_gene_pattern:

  Regular expression identifying non-marker gene mentions such as tRNA,
  ND, or ATP genes.

- raw:

  Raw data frame as read from a source file.

- source_name:

  Short label for the source, for example `"ncbi"` or `"bold"`.

- column_map:

  A named character vector where the names are the desired standardized
  column names and the values are the corresponding column names in
  `raw` (e.g., `c(species = "scientific_name")`).

- marker_fn:

  Function used to derive marker type from a description column. Default
  is `rb_extract_markers()`.

- description_col:

  Optional name of a free-text description column in `raw`.

- min_length:

  Minimum sequence length to retain. Default `1`.

## Value

- `rb_default_marker_patterns()` returns a named list of regular
  expressions.

- `rb_extract_markers()` returns a character string describing the
  detected marker and completeness, for example `"COI_partial"`,
  `"12S_unknown"`, `"genome_complete"`, `"multi_marker"`, or
  `"other_unknown"`.

- `rb_parse_source_table()` returns a standardized data frame.

## Details

The standardized schema used by barcurateR expects at least:

- `species`

- `sequence`

- `seq_type`

- `source`

Additional columns such as `sequence_id`, `description`, `length`, and
source-specific metadata may also be retained.

`rb_parse_source_table()` uses a user-supplied `column_map` to rename
raw columns into the standardized schema. If a description column is
provided, marker type is inferred using `rb_extract_markers()` and then
simplified into standardized `seq_type` labels (e.g., "`COI`",
"`other_unknown`").

## See also

Other parsing:
[`rb_ambiguity_utils`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md),
[`rb_combine_resolve`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md),
[`rb_sequence_prep`](https://zkua01.github.io/barcurateR/reference/rb_sequence_prep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
raw_ncbi <- read.csv("ncbi_sequences.csv", stringsAsFactors = FALSE)

parsed_ncbi <- rb_parse_source_table(
  raw_ncbi,
  source_name = "ncbi",
  column_map = c(
    sequence_id = "sequence_id",
    species = "species_query",
    sequence = "sequence"
  ),
  description_col = "description"
)
} # }
```
