# Write curated reference data to a SQLite database

Functions to write curated reference sequences, QC results, barcode gap
metrics, and ambiguity records into a SQLite database. These functions
do **not** filter data — downstream filtering is handled by
[`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md),
[`rb_filter_reference()`](https://zkua01.github.io/barcurateR/reference/rb_query.md),
and export functions.

- `rb_build_reference_db()`: One-step builder that creates a new SQLite
  file and writes all supplied tables in a single call.

- `rb_write_reference_table()`: Writes the main reference sequence
  table.

- `rb_write_qc_table()`: Writes the QC flag table.

- `rb_write_barcode_gap_table()`: Writes barcode gap metrics.

- `rb_write_ambiguity_table()`: Writes the pre-QC ambiguity archive.

## Usage

``` r
rb_write_reference_table(
  con,
  data,
  table_name = "reference_final",
  overwrite = TRUE,
  require_taxonomy = TRUE
)

rb_write_qc_table(con, data, table_name = "qc_reference", overwrite = TRUE)

rb_write_barcode_gap_table(
  con,
  data,
  table_name = "barcode_gap_metrics",
  overwrite = TRUE
)

rb_write_ambiguity_table(
  con,
  data,
  table_name = "ambiguous_sequences",
  overwrite = TRUE
)

rb_build_reference_db(
  path,
  reference_data,
  qc_data = NULL,
  barcode_gap_data = NULL,
  ambiguity_data = NULL,
  reference_table = "reference_final",
  qc_table = "qc_reference",
  barcode_gap_table = "barcode_gap_metrics",
  ambiguity_table = "ambiguous_sequences",
  overwrite = TRUE,
  require_taxonomy = TRUE
)
```

## Arguments

- con:

  A `DBIConnection` object to an open SQLite database.

- data:

  A data frame to write to the database.

- table_name:

  Name of the table to write.

- overwrite:

  Logical. Overwrite existing tables? Default `TRUE`.

- require_taxonomy:

  Logical. If `TRUE`, enforces that the reference table contains all
  seven taxonomic ranks (kingdom through species) before writing.
  Default `TRUE`.

- path:

  Path to the SQLite database file to create or overwrite.

- reference_data:

  Data frame of final curated reference sequences.

- qc_data:

  Optional data frame of QC results.

- barcode_gap_data:

  Optional data frame of barcode gap metrics.

- ambiguity_data:

  Optional data frame of ambiguity records.

- reference_table:

  Name for the reference table (default `"reference_final"`).

- qc_table:

  Name for the QC table (default `"qc_reference"`).

- barcode_gap_table:

  Name for the barcode gap table (default `"barcode_gap_metrics"`).

- ambiguity_table:

  Name for the ambiguity table (default `"ambiguous_sequences"`).

## Value

- `rb_build_reference_db()` returns the normalized path to the created
  database file (invisibly).

- The individual `rb_write_*` functions return `TRUE` invisibly on
  success.

## Details

The default table names are:

|                     |                       |
|---------------------|-----------------------|
| Table               | Default name          |
| Reference sequences | `reference_final`     |
| QC results          | `qc_reference`        |
| Barcode gap metrics | `barcode_gap_metrics` |
| Ambiguity archive   | `ambiguous_sequences` |

Before writing, the function automatically ensures that the following
columns exist (adding defaults where missing):

- `qc_flag` (default `"pass"`)

- `source` (default `"unknown"`)

- `seq_type` (standardized from `gene` or `marker` columns)

- `unique_code` (auto-generated from deduplicated sequences)

List-columns are rejected and factor columns are converted to character
before writing.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build a complete database in one call
rb_build_reference_db(
  path = "my_reference.db",
  reference_data = final_data,
  qc_data = qc_results,
  barcode_gap_data = gap_metrics
)

# Or write tables individually to an existing connection
con <- DBI::dbConnect(RSQLite::SQLite(), "my_reference.db")
rb_write_reference_table(con, final_data)
rb_write_qc_table(con, qc_results)
DBI::dbDisconnect(con)
} # }
```
