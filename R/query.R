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

rb_get_sequences <- function(con, species = NULL, genus = NULL, family = NULL,
                             marker = NULL, source = NULL, occurrence = NULL,
                             qc_flag = "pass", exact = FALSE) {
  conditions <- rb_sql_conditions(
    con, species = species, genus = genus, family = family,
    marker = marker, source = source, occurrence = occurrence, qc_flag = qc_flag, exact = exact
  )
  sql <- "select * from yzfishdb_final"
  if (length(conditions) > 0) {
    sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  }
  DBI::dbGetQuery(con, sql)
}

rb_list_taxa <- function(con, rank = "species", occurrence = NULL, marker = NULL, qc_flag = "pass", exact = FALSE) {
  valid <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (!rank %in% valid) stop("Unsupported rank: ", rank, call. = FALSE)
  rank_sql <- as.character(DBI::dbQuoteIdentifier(con, rank))
  conditions <- rb_sql_conditions(con, occurrence = occurrence, marker = marker, qc_flag = qc_flag, exact = exact)
  sql <- paste0("select ", rank_sql, " as ", rank, ", count(*) as n_sequences from yzfishdb_final")
  if (length(conditions) > 0) {
    sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  }
  sql <- paste(sql, "group by", rank_sql, "order by", rank_sql)
  DBI::dbGetQuery(con, sql)
}

#' Filter a reference data.frame by marker/source/occurrence/length/qc_flag
rb_filter_reference <- function(data, marker = NULL, source = NULL, occurrence = NULL,
                                 min_length = NULL, max_length = NULL,
                                 qc_flag = "pass", exact = FALSE) {
  out <- data
  if (!is.null(marker)) out <- out[out$seq_type %in% marker, , drop = FALSE]
  if (!is.null(source)) out <- out[out$source %in% source, , drop = FALSE]
  if (!is.null(occurrence)) out <- out[out$occurrence %in% occurrence, , drop = FALSE]
  if (!is.null(qc_flag)) {
    qc_flag <- qc_flag[!is.na(qc_flag)]
    keep <- rep(FALSE, nrow(out))
    for (flag in qc_flag) {
      if (exact) {
        keep <- keep | (out$qc_flag == flag)
      } else {
        if (identical(flag, "pass")) {
          keep <- keep | (out$qc_flag == "pass")
        } else {
          # BUGFIX 2026-08-03: was `out$qcflag` (nonexistent column —
          # missing underscore) instead of `out$qc_flag`. A nonexistent
          # column returns NULL, and grepl(flag, NULL) silently returns
          # logical(0), so `keep <- keep | logical(0)` left `keep`
          # unchanged with no error thrown. Net effect: every non-exact,
          # non-"pass" qc_flag exclusion silently did nothing in every
          # prior call. Fixed to reference the real column, out$qc_flag.
          keep <- keep | grepl(flag, out$qc_flag, fixed = TRUE)
        }
      }
    }
    out <- out[keep, , drop = FALSE]
  }

  seq_len <- nchar(rb_clean_sequence(out$sequence))
  if (!is.null(min_length)) out <- out[seq_len >= min_length, , drop = FALSE]
  if (!is.null(max_length)) out <- out[seq_len <= max_length, , drop = FALSE]
  row.names(out) <- NULL
  out
}

#' Exclude sequences matching any of a set of QC flags (in-memory filter)
#'
#' Thin wrapper around rb_filter_reference() for the common "keep
#' everything except these flags" case.
rb_qc_filter <- function(data, qc_flag_col = "qc_flag",
                         exclude = c("contaminant", "numt", "divergent",
                                     "frameshifted", "internal_stop",
                                     "sequence_short", "sequence_gap")) {
  if (!identical(qc_flag_col, "qc_flag")) {
    names(data)[names(data) == qc_flag_col] <- "qc_flag"
  }
  rb_filter_reference(data, qc_flag = exclude, exact = FALSE)
}


rb_taxonomy <- function(con, species = NULL, qc_flag = "pass", exact = FALSE) {
  conditions <- rb_sql_conditions(con, species = species, qc_flag = qc_flag, exact = exact)
  sql <- paste(
    "select distinct kingdom, phylum, class, `order`, family, genus, species, occurrence, habitat",
    "from yzfishdb_final"
  )
  if (length(conditions) > 0) sql <- paste(sql, "where", paste(conditions, collapse = " and "))
  DBI::dbGetQuery(con, sql)
}
