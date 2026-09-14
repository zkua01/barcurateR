#' Audit and inspect a curated reference database
#'
#' @description
#' Functions to generate diagnostic summaries and inspect the contents
#' of a curated reference database. These functions help users assess
#' data quality, marker coverage, source representation, and potential
#' issues before using the database for taxonomic assignment.
#'
#' **Database query functions** (require an open SQLite connection):
#'
#' * `rb_marker_coverage()`: Summarises sequences per taxonomic rank and marker.
#' * `rb_source_coverage()`: Summarises sequences per source and marker.
#' * `rb_qc_summary()`: Summarises QC flag distributions.
#' * `rb_barcode_gap()`: Retrieves barcode gap metrics from the database.
#' * `rb_ambiguity()`: Retrieves archived ambiguity records from the database.
#'
#' **File reader functions** (read from CSV files):
#'
#' * `rb_read_barcode_gap()`: Reads barcode gap metrics from a CSV file.
#' * `rb_read_ambiguity()`: Reads ambiguity records from a CSV file.
#'
#' @param con A `DBIConnection` object to an open barcurateR SQLite database.
#' @param rank Taxonomic rank to group by: `"species"`, `"genus"`, `"family"`,
#'   or `"order"`. Default `"species"`.
#' @param table_name Name of the reference table to query (default `"reference_final"`).
#' @param qc_table_name Name of the QC table to query (default `"qc_reference"`).
#' @param species Optional character vector of species names to filter by.
#' @param marker Optional character vector of marker names to filter by.
#' @param path Path to a CSV file. If `NULL`, a default bundled file is used.
#' @param action Optional filter for ambiguity resolution actions.
#'
#' @return
#' * `rb_marker_coverage()` returns a data frame with columns for the requested
#'   rank, `seq_type`, `n_sequences`, and `n_sources`.
#' * `rb_source_coverage()` returns a data frame with `source`, `seq_type`,
#'   `n_sequences`, and `n_species`.
#' * `rb_qc_summary()` returns a data frame with `qc_flag` and `n_sequences`,
#'   ordered by count descending.
#' * `rb_barcode_gap()` and `rb_read_barcode_gap()` return a data frame of
#'   barcode gap metrics.
#' * `rb_ambiguity()` and `rb_read_ambiguity()` return a data frame of
#'   ambiguity records.
#'
#' @details
#' All database query functions filter on `qc_flag = 'pass'` for coverage
#' functions (`rb_marker_coverage`, `rb_source_coverage`), ensuring only
#' passing sequences are counted. The `rb_qc_summary()` function reports
#' all QC flags without filtering.
#'
#' If the requested table does not exist in the database, `rb_barcode_gap()`
#' and `rb_ambiguity()` return an empty data frame rather than throwing an
#' error, allowing graceful degradation when optional tables are absent.
#'
#' @examples
#' \dontrun{
#' con <- rb_connect()
#'
#' # Marker coverage by species
#' rb_marker_coverage(con, rank = "species")
#'
#' # Source coverage
#' rb_source_coverage(con)
#'
#' # QC summary
#' rb_qc_summary(con)
#'
#' # Barcode gap metrics for a specific marker
#' rb_barcode_gap(con, marker = "COI")
#'
#' # Ambiguity records
#' rb_ambiguity(con)
#'
#' rb_disconnect(con)
#' }
#'
#' @name rb_diagnostics
#' @family diagnostics
NULL

#' @rdname rb_diagnostics
#' @export
rb_marker_coverage <- function(con, rank = "species", table_name = "reference_final") {
  valid <- c("species", "genus", "family", "order")
  if (!rank %in% valid) stop("Unsupported rank: ", rank, call. = FALSE)
  rank_sql <- as.character(DBI::dbQuoteIdentifier(con, rank))
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  sql <- paste0(
    "select ", rank_sql, " as ", rank, ", seq_type, count(*) as n_sequences, ",
    "count(distinct source) as n_sources from ", tbl,
    " where qc_flag = 'pass' group by ", rank_sql, ", seq_type ",
    "order by ", rank_sql, ", seq_type"
  )
  DBI::dbGetQuery(con, sql)
}

#' @rdname rb_diagnostics
#' @export
rb_source_coverage <- function(con, table_name = "reference_final") {
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  DBI::dbGetQuery(con, paste(
    "select source, seq_type, count(*) as n_sequences,",
    "count(distinct species) as n_species from", tbl,
    "where qc_flag = 'pass'",
    "group by source, seq_type order by source, seq_type"
  ))
}

#' @rdname rb_diagnostics
#' @export
rb_qc_summary <- function(con, table_name = "reference_final", qc_table_name = "qc_reference") {
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  qc_tbl <- DBI::dbQuoteIdentifier(con, qc_table_name)
  if (qc_table_name %in% DBI::dbListTables(con)) {
    DBI::dbGetQuery(con, paste(
      "select qc_flag, count(*) as n_sequences from", qc_tbl,
      "group by qc_flag order by n_sequences desc"
    ))
  } else {
    DBI::dbGetQuery(con, paste(
      "select qc_flag, count(*) as n_sequences from", tbl,
      "group by qc_flag order by n_sequences desc"
    ))
  }
}

#' @rdname rb_diagnostics
#' @export
rb_barcode_gap <- function(con, species = NULL, marker = NULL, table_name = "barcode_gap_metrics") {
  if (!table_name %in% DBI::dbListTables(con)) {
    return(data.frame())
  }
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  conditions <- character()
  if (!is.null(species)) conditions <- c(conditions, paste0("species in (", rb_quote_in(con, species), ")"))
  if (!is.null(marker)) conditions <- c(conditions, paste0("marker in (", rb_quote_in(con, marker), ")"))
  sql <- paste("select * from", tbl)
  if (length(conditions) > 0) sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  DBI::dbGetQuery(con, sql)
}

#' @rdname rb_diagnostics
#' @export
rb_ambiguity <- function(con, table_name = "ambiguous_sequences") {
  if (!table_name %in% DBI::dbListTables(con)) {
    return(data.frame())
  }
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  DBI::dbGetQuery(con, paste("select * from", tbl))
}

#' @rdname rb_diagnostics
#' @export
rb_read_barcode_gap <- function(path = NULL, species = NULL, marker = NULL) {
  if (is.null(path)) {
    path <- rb_default_data_file("barcode_gap_metrics_enhanced.csv")
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Barcode gap metrics file does not exist: ", path, call. = FALSE)
  }
  gap <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(species)) gap <- gap[gap$species %in% species, , drop = FALSE]
  if (!is.null(marker)) gap <- gap[gap$marker %in% marker, , drop = FALSE]
  row.names(gap) <- NULL
  gap
}

#' @rdname rb_diagnostics
#' @export
rb_read_ambiguity <- function(path = NULL, action = NULL) {
  if (is.null(path)) {
    path <- rb_default_data_file("ambiguous_sequences_tbl_cleaned.csv")
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Ambiguous sequence table does not exist: ", path, call. = FALSE)
  }
  amb <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(action)) amb <- amb[amb$action %in% action, , drop = FALSE]
  row.names(amb) <- NULL
  amb
}
