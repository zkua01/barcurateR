#' Query and filter a curated reference database
#'
#' @description
#' Functions to retrieve, filter, and inspect reference sequences from a
#' barcurateR SQLite database or an in-memory data frame.
#'
#' * `rb_get_sequences()`: Retrieves sequences from a SQLite database with
#'   optional taxonomic, marker, source, occurrence, and QC filters.
#' * `rb_list_taxa()`: Lists taxa at a given taxonomic rank with sequence
#'   counts.
#' * `rb_filter_reference()`: Filters an in-memory reference data frame
#'   using the same criteria as the database query functions.
#' * `rb_taxonomy()`: Retrieves distinct taxonomy records from the database.
#'
#' @param con A `DBIConnection` object to an open barcurateR SQLite database.
#' @param data A data frame containing reference sequences (for
#'   `rb_filter_reference()`).
#' @param species,genus,family Optional character vectors to filter by
#'   taxonomic rank.
#' @param marker Optional character vector of marker types to filter by
#'   (e.g., `"COI"`, `"12S"`). Maps to the `seq_type` column.
#' @param source Optional character vector of source names to filter by
#'   (e.g., `"ncbi"`, `"bold"`).
#' @param occurrence Optional character vector of occurrence categories
#'   to filter by (e.g., `"native"`, `"introduced"`).
#' @param qc_flag QC flag filter. Default `"pass"`. When `exact = FALSE`,
#'   non-`"pass"` values are matched using SQL `LIKE '%value%'`. When
#'   `exact = TRUE`, exact equality is used. Set to `NULL` to skip
#'   QC filtering.
#' @param exact Logical. If `TRUE`, use exact matching for `qc_flag`
#'   filters. If `FALSE` (default), use pattern matching for non-`"pass"`
#'   flags.
#' @param rank Taxonomic rank for `rb_list_taxa()`. One of `"kingdom"`,
#'   `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"`, or
#'   `"species"`. Default `"species"`.
#' @param table_name Name of the reference table to query.
#'   Default `"reference_final"`.
#' @param min_length,max_length Optional minimum and maximum sequence
#'   length filters for `rb_filter_reference()`.
#'
#' @return
#' * `rb_get_sequences()` returns a data frame of reference sequences
#'   matching the filter criteria.
#' * `rb_list_taxa()` returns a data frame with columns for the requested
#'   rank and `n_sequences`.
#' * `rb_filter_reference()` returns a filtered data frame.
#' * `rb_taxonomy()` returns a data frame of distinct taxonomy records
#'   with all seven ranks plus optional `occurrence` and `habitat` columns.
#'
#' @details
#' All database query functions use the generic table name
#' `"reference_final"` by default. This can be overridden via the
#' `table_name` parameter.
#'
#' `rb_get_sequences()` and `rb_list_taxa()` filter on `qc_flag = 'pass'`
#' by default, ensuring only quality-checked sequences are returned.
#' Set `qc_flag = NULL` to retrieve all sequences regardless of QC status.
#'
#' `rb_taxonomy()` requires all seven taxonomic ranks (`kingdom` through
#' `species`) to be present in the table. The columns `occurrence` and
#' `habitat` are included in the output if they exist in the table.
#'
#' @examples
#' \dontrun{
#' con <- rb_connect()
#'
#' # Get all passing COI sequences
#' coi <- rb_get_sequences(con, marker = "COI")
#'
#' # List all species with sequence counts
#' taxa <- rb_list_taxa(con, rank = "species")
#'
#' # Get taxonomy for a specific species
#' tax <- rb_taxonomy(con, species = "Danio rerio")
#'
#' # Filter an in-memory data frame
#' filtered <- rb_filter_reference(refs, marker = "12S", min_length = 100)
#'
#' rb_disconnect(con)
#' }
#'
#' @name rb_query
#' @family database queries
NULL

#' Build SQL WHERE conditions for reference queries
#'
#' Internal helper that constructs SQL filter clauses from named
#' parameters. Used by `rb_get_sequences()`, `rb_list_taxa()`, and
#' `rb_taxonomy()`.
#'
#' @keywords internal
#' @noRd
rb_sql_conditions <- function(con, species = NULL, genus = NULL, family = NULL,
                              marker = NULL, source = NULL, occurrence = NULL,
                              qc_flag = NULL, exact = FALSE) {
  conditions <- character()
  add_in <- function(field, values) {
    if (is.null(values)) return(NULL)
    quoted_field <- as.character(DBI::dbQuoteIdentifier(con, field))
    paste0(quoted_field, " in (", rb_quote_in(con, values), ")")
  }
  conditions <- c(
    conditions,
    add_in("species", species),
    add_in("genus", genus),
    add_in("family", family),
    add_in("seq_type", marker),
    add_in("source", source),
    add_in("occurrence", occurrence)
  )
  if (!is.null(qc_flag)) {
    qc_flag <- qc_flag[!is.na(qc_flag)] # Remove NAs just in case
    qc_conds <- character()
    
    for (flag in qc_flag) {
      if(exact) {
        quoted_flag <- DBI::dbQuoteString(con, flag)
        qc_conds <- c(qc_conds, paste0("qc_flag = ", quoted_flag))
      } else{
        if (identical(flag, "pass")) {
          qc_conds <- c(qc_conds, "qc_flag = 'pass'")
        } else {
          pattern <- DBI::dbQuoteString(con, paste0("%", flag, "%"))
          qc_conds <- c(qc_conds, paste0("qc_flag like ", pattern))
        }
      }
    }
    
    # If multiple flags are provided, join them with OR and wrap in parentheses
    if (length(qc_conds) > 1) {
      conditions <- c(conditions, paste0("(", paste(qc_conds, collapse = " or "), ")"))
    } else if (length(qc_conds) == 1) {
      conditions <- c(conditions, qc_conds)
    }
  }
  conditions[!is.na(conditions) & nzchar(conditions)]
}

#' @rdname rb_query
#' @export
rb_get_sequences <- function(con, species = NULL, genus = NULL, family = NULL,
                             marker = NULL, source = NULL, occurrence = NULL,
                             qc_flag = "pass", exact = FALSE,
                             table_name = "reference_final") {
  conditions <- rb_sql_conditions(
    con, species = species, genus = genus, family = family,
    marker = marker, source = source, occurrence = occurrence, qc_flag = qc_flag, exact = exact
  )
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  sql <- paste("select * from", tbl)
  if (length(conditions) > 0) {
    sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  }
  DBI::dbGetQuery(con, sql)
}

#' @rdname rb_query
#' @export
rb_list_taxa <- function(con, rank = "species", occurrence = NULL, marker = NULL, qc_flag = "pass", exact = FALSE, table_name = "reference_final") {
  valid <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (!rank %in% valid) stop("Unsupported rank: ", rank, call. = FALSE)
  rank_sql <- DBI::dbQuoteIdentifier(con, rank)
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  conditions <- rb_sql_conditions(con, occurrence = occurrence, marker = marker, qc_flag = qc_flag, exact = exact)
  sql <- paste0("select ", rank_sql, " as ", rank, ", count(*) as n_sequences from ", tbl)
  if (length(conditions) > 0) {
    sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  }
  sql <- paste(sql, "group by", rank_sql, "order by", rank_sql)
  DBI::dbGetQuery(con, sql)
}

#' @rdname rb_query
#' @export
rb_filter_reference <- function(data, marker = NULL, source = NULL, occurrence = NULL, min_length = NULL, max_length = NULL, qc_flag = "pass", exact = FALSE) {
  out <- data
  if (!is.null(marker)) out <- out[out$seq_type %in% marker, , drop = FALSE]
  if (!is.null(source)) out <- out[out$source %in% source, , drop = FALSE]
  if (!is.null(occurrence)) out <- out[out$occurrence %in% occurrence, , drop = FALSE]
  if (!is.null(qc_flag)){
    qc_flag <- qc_flag[!is.na(qc_flag)]
    keep <- rep(FALSE, nrow(out))
    for(flag in qc_flag){
      if(exact){
        keep <- keep | (out$qc_flag == flag)
      } else{
        if(identical(flag, "pass")){
          keep <- keep | (out$qc_flag == "pass")
        } else {
          keep <- keep | grepl(flag, out$qc_flag, fixed = TRUE)
        }
      }
    }
    out <- out[keep,,drop=FALSE]
  }
  
  seq_len <- nchar(rb_clean_sequence(out$sequence))
  if (!is.null(min_length)) out <- out[seq_len >= min_length, , drop = FALSE]
  if (!is.null(max_length)) out <- out[seq_len <= max_length, , drop = FALSE]
  row.names(out) <- NULL
  out
}

#' @rdname rb_query
#' @export
rb_taxonomy <- function(con, species = NULL, qc_flag = "pass", exact = FALSE, table_name = "reference_final") {
  if (!table_name %in% DBI::dbListTables(con)) {
    stop("Table not found: ", table_name, call. = FALSE)
  }
  fields <- DBI::dbListFields(con, table_name)
  tax_ranks <- c("kingdom","phylum","class","order","family","genus","species")
  missing_tax <- setdiff(tax_ranks, fields)
  if (length(missing_tax)>0){
    stop("Reference table is missing required taxonomy columns: ",
         paste(missing_tax, collapse = ", "),
         call. = FALSE)
  }
  optional_cols <- intersect(c("occurrence","habitat"), fields)
  select_cols <- c(tax_ranks, optional_cols)
  
  use_qc <- !is.null(qc_flag) && "qc_flag" %in% fields
  
  conditions <- rb_sql_conditions(con, species = species, qc_flag = if (use_qc) qc_flag else NULL, exact = exact)
  tbl <- DBI::dbQuoteIdentifier(con, table_name)
  order_col <- DBI::dbQuoteIdentifier(con, "order")
  sql <- paste(
    "select distinct", paste(as.character(DBI::dbQuoteIdentifier(con, select_cols)), collapse = ", "), "from", tbl
  )
  if (length(conditions) > 0) sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  DBI::dbGetQuery(con, sql)
}
