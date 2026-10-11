# Query and filter a curated reference database

Functions to retrieve, filter, and inspect reference sequences from a
barcurateR SQLite database or an in-memory data frame.

- `rb_get_sequences()`: Retrieves sequences from a SQLite database with
  optional taxonomic, marker, source, occurrence, and QC filters.

- `rb_list_taxa()`: Lists taxa at a given taxonomic rank with sequence
  counts.

- `rb_filter_reference()`: Filters an in-memory reference data frame
  using the same criteria as the database query functions.

- `rb_taxonomy()`: Retrieves distinct taxonomy records from the
  database.

## Usage

``` r
rb_get_sequences(
  con,
  species = NULL,
  genus = NULL,
  family = NULL,
  marker = NULL,
  source = NULL,
  occurrence = NULL,
  qc_flag = "pass",
  exact = FALSE,
  table_name = "reference_final"
)

rb_list_taxa(
  con,
  rank = "species",
  occurrence = NULL,
  marker = NULL,
  qc_flag = "pass",
  exact = FALSE,
  table_name = "reference_final"
)

rb_filter_reference(
  data,
  marker = NULL,
  source = NULL,
  occurrence = NULL,
  min_length = NULL,
  max_length = NULL,
  qc_flag = "pass",
  exact = FALSE
)

rb_taxonomy(
  con,
  species = NULL,
  qc_flag = "pass",
  exact = FALSE,
  table_name = "reference_final"
)
```

## Arguments

- con:

  A `DBIConnection` object to an open barcurateR SQLite database.

- species, genus, family:

  Optional character vectors to filter by taxonomic rank.

- marker:

  Optional character vector of marker types to filter by (e.g., `"COI"`,
  `"12S"`). Maps to the `seq_type` column.

- source:

  Optional character vector of source names to filter by (e.g.,
  `"ncbi"`, `"bold"`).

- occurrence:

  Optional character vector of occurrence categories to filter by (e.g.,
  `"native"`, `"introduced"`). Requires `occurrence` column to be
  present.

- qc_flag:

  QC flag filter. Default `"pass"`. When `exact = FALSE`, non-`"pass"`
  values are matched using SQL `LIKE '%value%'`. When `exact = TRUE`,
  exact equality is used. Set to `NULL` to skip QC filtering.

- exact:

  Logical. If `TRUE`, use exact matching for `qc_flag` filters. If
  `FALSE` (default), use pattern matching for non-`"pass"` flags.

- table_name:

  Name of the reference table to query. Default `"reference_final"`.

- rank:

  Taxonomic rank for `rb_list_taxa()`. One of `"kingdom"`, `"phylum"`,
  `"class"`, `"order"`, `"family"`, `"genus"`, or `"species"`. Default
  `"species"`.

- data:

  A data frame containing reference sequences (for
  `rb_filter_reference()`).

- min_length, max_length:

  Optional minimum and maximum sequence length filters for
  `rb_filter_reference()`.

## Value

- `rb_get_sequences()` returns a data frame of reference sequences
  matching the filter criteria.

- `rb_list_taxa()` returns a data frame with columns for the requested
  rank and `n_sequences`.

- `rb_filter_reference()` returns a filtered data frame.

- `rb_taxonomy()` returns a data frame of distinct taxonomy records with
  all seven ranks plus optional `occurrence` and `habitat` columns.

## Details

All database query functions use the generic table name
`"reference_final"` by default. This can be overridden via the
`table_name` parameter.

`rb_get_sequences()` and `rb_list_taxa()` filter on `qc_flag = 'pass'`
by default, ensuring only quality-checked sequences are returned. Set
`qc_flag = NULL` to retrieve all sequences regardless of QC status.

`rb_taxonomy()` requires all seven taxonomic ranks (`kingdom` through
`species`) to be present in the table. The columns `occurrence` and
`habitat` are included in the output if they exist in the table.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- rb_connect()

# Get all passing COI sequences
coi <- rb_get_sequences(con, marker = "COI")

# List all species with sequence counts
taxa <- rb_list_taxa(con, rank = "species")

# Get taxonomy for a specific species
tax <- rb_taxonomy(con, species = "Danio rerio")

# Filter an in-memory data frame
filtered <- rb_filter_reference(refs, marker = "12S", min_length = 100)

rb_disconnect(con)
} # }
```
