rb_marker_coverage <- function(con, rank = "species", table_name = "yzfishdb_final") {
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

rb_source_coverage <- function(con, table_name = "yzfishdb_final") {
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  DBI::dbGetQuery(con, paste(
    "select source, seq_type, count(*) as n_sequences,",
    "count(distinct species) as n_species from", tbl,
    "where qc_flag = 'pass'",
    "group by source, seq_type order by source, seq_type"
  ))
}

rb_qc_summary <- function(con, table_name = "yzfishdb_final", qc_table_name = "qc_reference_p2") {
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

rb_ambiguity <- function(con, table_name = "ambiguous_sequences") {
  if (!table_name %in% DBI::dbListTables(con)) {
    return(data.frame())
  }
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  DBI::dbGetQuery(con, paste("select * from", tbl))
}

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
