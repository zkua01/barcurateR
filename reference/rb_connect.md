# Connect to and manage barcurateR SQLite databases

Functions to connect to, inspect, and manage barcurateR SQLite reference
databases. These functions provide the primary interface for interacting
with curated reference sequence databases.

- `rb_connect()`: Opens a connection to a SQLite database.

- `rb_disconnect()`: Safely closes an active database connection.

- `rb_tables()`: Lists all tables in the database with row counts.

- `rb_read_table()`: Reads an entire table into a data frame.

- `rb_available_datasets()`: Lists all available bundled, cached, or
  legacy databases.

## Usage

``` r
rb_available_datasets()

rb_connect(path = NULL, download = FALSE, url = rb_db_url())

rb_disconnect(con)

rb_tables(con)

rb_read_table(con, table)

rb_available_datasets()
```

## Arguments

- path:

  Path to a SQLite database file. If `NULL`, the default configured
  database is resolved (downloaded if requested, or the first available
  local database is used).

- download:

  Logical. If `TRUE` and the database file does not exist, attempt to
  download it from the default URL.

- url:

  The URL to download the database from (used when `download = TRUE`).

- con:

  A `DBIConnection` object returned by `rb_connect()`.

- table:

  The name of the table to read.

## Value

- `rb_connect()` returns a `DBIConnection` object (S4 class
  `SQLiteConnection`).

- `rb_disconnect()` returns `TRUE` invisibly.

- `rb_tables()` returns a data frame with columns `table` and `n_rows`.

- `rb_read_table()` returns a data frame containing the table contents.

- `rb_available_datasets()` returns a data frame with columns `name`,
  `version`, and `path`. It searches the package's `extdata` directory,
  the `barcurateR` cache directory, and the legacy `regionbarcoder`
  cache directory for any `.db` or `.sqlite` files.

## Details

The default database resolution order in `rb_connect()` when
`path = NULL` is:

1.  If `download = TRUE` and no cached file exists, download from
    [`rb_db_url()`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md).

2.  Otherwise, use the best available local database (preferring
    `full-local` YZFishDB over `demo` or `user-curated` databases).

## See also

Other database management:
[`rb_db_install`](https://zkua01.github.io/barcurateR/reference/rb_db_install.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# List all available local databases
avail <- rb_available_datasets()

# Connect to the default database
con <- rb_connect()
rb_tables(con)

# Read a specific table
refs <- rb_read_table(con, "reference_final")

# Close the connection
rb_disconnect(con)
} # }
```
