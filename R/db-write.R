#' Write curated reference data to a SQLite database
#'
#' @description
#' Functions to write curated reference sequences, QC results, barcode gap
#' metrics, and ambiguity records into a SQLite database. These functions
#' do **not** filter data — downstream filtering is handled by
#' [rb_get_sequences()], [rb_filter_reference()], and export functions.
#'
#' * `rb_build_reference_db()`: One-step builder that creates a new SQLite
#'   file and writes all supplied tables in a single call.
#' * `rb_write_reference_table()`: Writes the main reference sequence table.
#' * `rb_write_qc_table()`: Writes the QC flag table.
#' * `rb_write_barcode_gap_table()`: Writes barcode gap metrics.
#' * `rb_write_ambiguity_table()`: Writes the pre-QC ambiguity archive.
#'
#' @param path Path to the SQLite database file to create or overwrite.
#' @param con A `DBIConnection` object to an open SQLite database.
#' @param data A data frame to write to the database.
#' @param reference_data Data frame of final curated reference sequences.
#' @param qc_data Optional data frame of QC results.
#' @param barcode_gap_data Optional data frame of barcode gap metrics.
#' @param ambiguity_data Optional data frame of ambiguity records.
#' @param table_name Name of the table to write.
#' @param reference_table Name for the reference table (default `"reference_final"`).
#' @param qc_table Name for the QC table (default `"qc_reference"`).
#' @param barcode_gap_table Name for the barcode gap table (default `"barcode_gap_metrics"`).
#' @param ambiguity_table Name for the ambiguity table (default `"ambiguous_sequences"`).
#' @param overwrite Logical. Overwrite existing tables? Default `TRUE`.
#' @param require_taxonomy Logical. If `TRUE`, enforces that the reference
#'   table contains all seven taxonomic ranks (kingdom through species)
#'   before writing. Default `TRUE`.
#'
#' @return
#' * `rb_build_reference_db()` returns the normalized path to the created
#'   database file (invisibly).
#' * The individual `rb_write_*` functions return `TRUE` invisibly on success.
#'
#' @details
#' The default table names are:
#'
#' | Table | Default name |
#' |---|---|
#' | Reference sequences | `reference_final` |
#' | QC results | `qc_reference` |
#' | Barcode gap metrics | `barcode_gap_metrics` |
#' | Ambiguity archive | `ambiguous_sequences` |
#'
#' Before writing, `rb_prepare_reference_table()` ensures that the following
#' columns exist (adding defaults where missing):
#'
#' * `qc_flag` (default `"pass"`)
#' * `source` (default `"unknown"`)
#' * `seq_type` (standardized from `gene` or `marker` columns)
#' * `unique_code` (auto-generated from deduplicated sequences)
#'
#' List-columns are rejected and factor columns are converted to character
#' before writing.
#'
#' @examples
#' \dontrun{
#' # Build a complete database in one call
#' rb_build_reference_db(
#'   path = "my_reference.db",
#'   reference_data = final_data,
#'   qc_data = qc_results,
#'   barcode_gap_data = gap_metrics
#' )
#'
#' # Or write tables individually to an existing connection
#' con <- DBI::dbConnect(RSQLite::SQLite(), "my_reference.db")
#' rb_write_reference_table(con, final_data)
#' rb_write_qc_table(con, qc_results)
#' DBI::dbDisconnect(con)
#' }
#'
#' @name rb_write_db
#' @family database writing
NULL

#' Internal: prepare a data frame for SQLite writing
#'
#' Converts factors to character and rejects list-columns.
#'
#' @keywords internal
#' @noRd
rb_prepare_for_sqlite <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  
  list_cols <- vapply(data, is.list, logical(1))
  
  if (any(list_cols)) {
    stop(
      "Cannot write list column(s) to SQLite: ",
      paste(names(data)[list_cols], collapse = ", "),
      call. = FALSE
    )
  }
  
  factor_cols <- vapply(data, is.factor, logical(1))
  data[factor_cols] <- lapply(data[factor_cols], as.character)
  
  data
}

#' Internal: ensure minimal reference table schema
#'
#' Adds missing `qc_flag`, `source`, `seq_type`, and `unique_code` columns.
#' Optionally validates taxonomy requirements.
#'
#' @keywords internal
#' @noRd
rb_prepare_reference_table <- function(data, require_taxonomy = TRUE) {
  rb_required_columns(data, c("species", "sequence"))
  
  # Ensure qc_flag exists so rb_get_sequences() can filter by it.
  if (!"qc_flag" %in% names(data)) {
    data$qc_flag <- "pass"
  }
  
  # Ensure source exists so rb_source_coverage() works.
  if (!"source" %in% names(data)) {
    data$source <- "unknown"
  }
  
  # Standardize marker column to seq_type.
  if (!"seq_type" %in% names(data)) {
    if ("gene" %in% names(data)) {
      data$seq_type <- data$gene
    } else if ("marker" %in% names(data)) {
      data$seq_type <- data$marker
    } else {
      data$seq_type <- "other_unknown"
    }
  }
  
  # Ensure unique_code exists for assignment functions.
  if (!"unique_code" %in% names(data)) {
    seq_clean <- rb_clean_sequence(data$sequence)
    unique_seqs <- unique(seq_clean)
    codes <- sprintf("seq_%05d", seq_along(unique_seqs))
    data$unique_code <- codes[match(seq_clean, unique_seqs)]
  }
  
  if (isTRUE(require_taxonomy)){
    rb_required_taxonomy_columns(data)
  }
  
  data
}

#' Internal: generic table writer
#'
#' Validates connection and data, then writes via [DBI::dbWriteTable()].
#'
#' @keywords internal
#' @noRd
rb_write_table <- function(con, table_name, data, overwrite = TRUE) {
  if (!DBI::dbIsValid(con)) {
    stop("Database connection is not valid.", call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  
  if (!nzchar(table_name)) {
    stop("`table_name` must be a non-empty string.", call. = FALSE)
  }
  
  data <- rb_prepare_for_sqlite(data)
  
  DBI::dbWriteTable(con, table_name, data, overwrite = overwrite)
  
  invisible(TRUE)
}

#' @rdname rb_write_db
#' @export
rb_write_reference_table <- function(con, data,
                                     table_name = "reference_final",
                                     overwrite = TRUE,
                                     require_taxonomy = TRUE) {
  data <- rb_prepare_reference_table(data, require_taxonomy = require_taxonomy)
  
  rb_write_table(
    con,
    table_name = table_name,
    data = data,
    overwrite = overwrite
  )
}

#' @rdname rb_write_db
#' @export
rb_write_qc_table <- function(con, data,
                              table_name = "qc_reference",
                              overwrite = TRUE) {
  
  if (!"qc_flag" %in% names(data)) {
    data$qc_flag <- "pass"
  }
  
  rb_write_table(
    con,
    table_name = table_name,
    data = data,
    overwrite = overwrite
  )
}

#' @rdname rb_write_db
#' @export
rb_write_barcode_gap_table <- function(con, data,
                                       table_name = "barcode_gap_metrics",
                                       overwrite = TRUE) {
  rb_required_columns(data, c("species", "marker"))
  
  rb_write_table(
    con,
    table_name = table_name,
    data = data,
    overwrite = overwrite
  )
}

#' @rdname rb_write_db
#' @export
rb_write_ambiguity_table <- function(con, data,
                                     table_name = "ambiguous_sequences",
                                     overwrite = TRUE) {
  rb_write_table(
    con,
    table_name = table_name,
    data = data,
    overwrite = overwrite
  )
}

#' @rdname rb_write_db
#' @export
rb_build_reference_db <- function(path,
                                  reference_data,
                                  qc_data = NULL,
                                  barcode_gap_data = NULL,
                                  ambiguity_data = NULL,
                                  reference_table = "reference_final",
                                  qc_table = "qc_reference",
                                  barcode_gap_table = "barcode_gap_metrics",
                                  ambiguity_table = "ambiguous_sequences",
                                  overwrite = TRUE,
                                  require_taxonomy = TRUE) {
  
  if (!nzchar(path)) {
    stop("`path` must be a non-empty file path.", call. = FALSE)
  }
  
  dest_dir <- dirname(path)
  
  if (nzchar(dest_dir) && !dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  rb_write_reference_table(
    con,
    data = reference_data,
    table_name = reference_table,
    overwrite = overwrite,
    require_taxonomy = require_taxonomy
  )
  
  if (!is.null(qc_data)) {
    rb_write_qc_table(
      con,
      data = qc_data,
      table_name = qc_table,
      overwrite = overwrite
    )
  }
  
  if (!is.null(barcode_gap_data)) {
    rb_write_barcode_gap_table(
      con,
      data = barcode_gap_data,
      table_name = barcode_gap_table,
      overwrite = overwrite
    )
  }
  
  if (!is.null(ambiguity_data)) {
    rb_write_ambiguity_table(
      con,
      data = ambiguity_data,
      table_name = ambiguity_table,
      overwrite = overwrite
    )
  }
  
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}