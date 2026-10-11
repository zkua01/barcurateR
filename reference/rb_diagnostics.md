# Audit and inspect a curated reference database

Functions to generate diagnostic summaries and inspect the contents of a
curated reference database. These functions help users assess data
quality, marker coverage, source representation, and potential issues
before using the database for taxonomic assignment.

**Database query functions** (require an open SQLite connection):

- `rb_marker_coverage()`: Summarises sequences per taxonomic rank and
  marker.

- `rb_source_coverage()`: Summarises sequences per source and marker.

- `rb_qc_summary()`: Summarises QC flag distributions.

- `rb_barcode_gap()`: Retrieves barcode gap metrics from the database.

- `rb_ambiguity()`: Retrieves archived ambiguity records from the
  database.

**File reader functions** (read from CSV files):

- `rb_read_barcode_gap()`: Reads barcode gap metrics from a CSV file.

- `rb_read_ambiguity()`: Reads ambiguity records from a CSV file.

## Usage

``` r
rb_marker_coverage(con, rank = "species", table_name = "reference_final")

rb_source_coverage(con, table_name = "reference_final")

rb_qc_summary(
  con,
  table_name = "reference_final",
  qc_table_name = "qc_reference"
)

rb_barcode_gap(
  con,
  species = NULL,
  marker = NULL,
  table_name = "barcode_gap_metrics"
)

rb_ambiguity(con, table_name = "ambiguous_sequences")

rb_read_barcode_gap(path = NULL, species = NULL, marker = NULL)

rb_read_ambiguity(path = NULL, action = NULL)
```

## Arguments

- con:

  A `DBIConnection` object to an open barcurateR SQLite database.

- rank:

  Taxonomic rank to group by: `"species"`, `"genus"`, `"family"`, or
  `"order"`. Default `"species"`.

- table_name:

  Name of the reference table to query (default `"reference_final"`).

- qc_table_name:

  Name of the QC table to query (default `"qc_reference"`).

- species:

  Optional character vector of species names to filter by.

- marker:

  Optional character vector of marker names to filter by.

- path:

  Path to a CSV file. If `NULL`, a default bundled file is used.

- action:

  Optional filter for ambiguity resolution actions.

## Value

- `rb_marker_coverage()` returns a data frame with columns for the
  requested rank, `seq_type`, `n_sequences`, and `n_sources`.

- `rb_source_coverage()` returns a data frame with `source`, `seq_type`,
  `n_sequences`, and `n_species`.

- `rb_qc_summary()` returns a data frame with `qc_flag` and
  `n_sequences`, ordered by count descending.

- `rb_barcode_gap()` and `rb_read_barcode_gap()` return a data frame of
  barcode gap metrics.

- `rb_ambiguity()` and `rb_read_ambiguity()` return a data frame of
  ambiguity records.

## Details

All database query functions filter on `qc_flag = 'pass'` for coverage
functions (`rb_marker_coverage`, `rb_source_coverage`), ensuring only
passing sequences are counted. The `rb_qc_summary()` function reports
all QC flags without filtering.

If the requested table does not exist in the database,
`rb_barcode_gap()` and `rb_ambiguity()` return an empty data frame
rather than throwing an error, allowing graceful degradation when
optional tables are absent.

If the separate QC table (`qc_table_name`) does not exist in the
database, `rb_qc_summary()` falls back to summarizing the `qc_flag`
column in the main reference table (`table_name`).

## Examples

``` r
if (FALSE) { # \dontrun{
con <- rb_connect()

# Marker coverage by species
rb_marker_coverage(con, rank = "species")

# Source coverage
rb_source_coverage(con)

# QC summary
rb_qc_summary(con)

# Barcode gap metrics for a specific marker
rb_barcode_gap(con, marker = "COI")

# Ambiguity records
rb_ambiguity(con)

rb_disconnect(con)
} # }
```
