#' Connect to and manage barcurateR SQLite databases
#'
#' @description
#' Functions to connect to, inspect, and manage barcurateR SQLite reference
#' databases. These functions provide the primary interface for interacting
#' with curated reference sequence databases.
#'
#' * `rb_connect()`: Opens a connection to a SQLite database.
#' * `rb_disconnect()`: Safely closes an active database connection.
#' * `rb_tables()`: Lists all tables in the database with row counts.
#' * `rb_read_table()`: Reads an entire table into a data frame.
#' * `rb_available_datasets()`: Lists available bundled or local databases.
#'
#' @param path Path to a SQLite database file. If `NULL`, the default
#'   bundled demo database is used.
#' @param download Logical. If `TRUE` and the database file does not exist,
#'   attempt to download it from the default URL.
#' @param url The URL to download the database from (used when `download = TRUE`).
#' @param con A `DBIConnection` object returned by `rb_connect()`.
#' @param table The name of the table to read.
#'
#' @return
#' * `rb_connect()` returns a `DBIConnection` object (S4 class `SQLiteConnection`).
#' * `rb_disconnect()` returns `TRUE` invisibly.
#' * `rb_tables()` returns a data frame with columns `table` and `n_rows`.
#' * `rb_read_table()` returns a data frame containing the table contents.
#' * `rb_available_datasets()` returns a data frame with columns `name`, `version`, and `path`.
#'
#' @details
#' The default database resolution order in `rb_connect()` when `path = NULL` is:
#' 1. If `download = TRUE` and no cached file exists, download from `rb_db_url()`.
#' 2. Otherwise, use the first path from `rb_available_datasets()`.
#'
#' @examples
#' \dontrun{
#' # Connect to the bundled demo database
#' con <- rb_connect()
#' rb_tables(con)
#'
#' # Read a specific table
#' refs <- rb_read_table(con, "reference_final")
#'
#' # Close the connection
#' rb_disconnect(con)
#' }
#'
#' @name rb_connect
#' @family database management
NULL

#' @rdname rb_connect
#' @export
rb_available_datasets <- function() {
  demo_path <- system.file("extdata", "small_refdb.sqlite", package = "barcurateR")
  if (!nzchar(demo_path)) {
    demo_path <- rb_find_full_yzfishdb()
  }
  if (!nzchar(demo_path)) {
    demo_path <- rb_demo_db_path()
  }
  data.frame(
    name = "small_refdb",
    version = if (grepl("YZFishDB[.]db$", demo_path)) "full-local" else "demo",
    path = demo_path,
    stringsAsFactors = FALSE
  )
}

#' @rdname rb_connect
#' @export
rb_connect <- function(path = NULL, download = FALSE, url = rb_db_url()) {
  if (is.null(path)) {
    if (isTRUE(download) && !file.exists(rb_db_path())) {
      path <- rb_install_db(url = url)
    } else {
      path <- rb_available_datasets()$path[[1]]
    }
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Database file does not exist: ", path, call. = FALSE)
  }
  DBI::dbConnect(RSQLite::SQLite(), normalizePath(path, winslash = "/", mustWork = TRUE))
}

#' @rdname rb_connect
#' @export
rb_disconnect <- function(con) {
  if (DBI::dbIsValid(con)) {
    DBI::dbDisconnect(con)
  }
  invisible(TRUE)
}

#' @rdname rb_connect
#' @export
rb_tables <- function(con) {
  tables <- DBI::dbListTables(con)
  counts <- vapply(tables, function(tab) {
    quoted <- as.character(DBI::dbQuoteIdentifier(con, tab))
    DBI::dbGetQuery(con, paste0("select count(*) as n from ", quoted))[[1]]
  }, numeric(1))
  data.frame(table = tables, n_rows = as.integer(counts), stringsAsFactors = FALSE)
}

#' @rdname rb_connect
#' @export
rb_read_table <- function(con, table) {
  available <- DBI::dbListTables(con)
  if (!table %in% available) {
    stop("Table not found: ", table, call. = FALSE)
  }
  DBI::dbReadTable(con, table)
}
