# ============================================================
# MODULE — SQLite database writing
#
# DESTINATION: R/db-write.R
#
# This module writes curated reference data, QC results, and
# barcode gap metrics into a SQLite database.
#
# It deliberately does not implement the original final
# filtering step. Downstream filtering can be handled by
# rb_get_sequences(), rb_filter_reference(), and export
# functions.
# ============================================================


# ------------------------------------------------------------
# Internal helper: prepare data frame for SQLite writing
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Internal helper: ensure minimal reference table schema
#
# This does not filter data. It only ensures that the columns
# needed by common downstream functions exist.
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Generic table writer
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Write main reference table
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Write QC table
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Write barcode gap metrics
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Optional: write ambiguity table
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# One-step database builder
# ------------------------------------------------------------

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