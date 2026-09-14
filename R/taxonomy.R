#' Taxonomy standardization, joining, and string formatting
#'
#' @description
#' Functions to standardize species names, join taxonomic lineage onto
#' reference data, and format taxonomy strings for export.
#'
#' * `rb_build_taxonomy_string()`: Formats the seven taxonomic ranks into
#'   a semicolon-delimited string, optionally with rank prefixes.
#' * `rb_standardize_species()`: Maps alternative species names to valid
#'   names using a reference table, with an optional external taxonomy
#'   lookup callback as fallback.
#' * `rb_join_taxonomy()`: Joins a taxonomy table onto reference data,
#'   filling missing genus/family/order with a configurable label and
#'   applying constant higher-rank defaults.
#'
#' @param data A data frame containing at least the seven taxonomic
#'   ranks (`kingdom` through `species`) for `rb_build_taxonomy_string()`,
#'   or at least the column specified by `by` for `rb_join_taxonomy()`.
#' @param style Formatting style for `rb_build_taxonomy_string()`.
#'   One of `"plain"` (default) or `"rank_prefix"`.
#' @param names Character vector of species names to standardize.
#' @param reference_table A data frame mapping alternative names to
#'   valid names for `rb_standardize_species()`.
#' @param taxonomy_table A data frame containing taxonomic lineage to join
#'   onto `data` for `rb_join_taxonomy()`. Must contain the column specified by `by`, and optionally
#'   columns for `genus`, `family`, `order`, and higher ranks.
#' @param name_col Column in `reference_table` holding alternative/known
#'   names. Default `"alternative_clean"`.
#' @param valid_col Column in `reference_table` holding the valid name.
#'   Default `"valid_name"`.
#' @param taxonomy_lookup Optional callback `function(name) -> valid_name`
#'   used when `reference_table` has no match (e.g., an external catalog
#'   API). `NULL` by default.
#' @param by Column name to join on. Default `"species"`.
#' @param defaults Named list of constant higher-rank values applied to
#'   every row, e.g., `list(kingdom = "Metazoa", phylum = "Chordata")`.
#' @param unassigned_label Label used when genus/family/order cannot be
#'   resolved. Default `"UNASSIGNED"`.
#'
#' @return
#' * `rb_build_taxonomy_string()` returns a character vector of formatted
#'   taxonomy strings, one per row.
#' * `rb_standardize_species()` returns a character vector of standardized
#'   species names, with `NA` for unmatched names.
#' * `rb_join_taxonomy()` returns the input data frame with taxonomic
#'   lineage columns joined.
#'
#' @details
#' `rb_build_taxonomy_string()` supports two styles:
#'
#' | Style | Example output |
#' |---|---|
#' | `"plain"` | `Animalia;Chordata;Actinopterygii;...` |
#' | `"rank_prefix"` | `k__Animalia;p__Chordata;c__Actinopterygii;...` |
#'
#' Missing or empty values are replaced with `"unassigned"`.
#'
#' `rb_join_taxonomy()` guards against silent Cartesian joins by checking
#' for duplicate keys in `taxonomy_table` before merging. If duplicates
#' are found, it stops with an informative error.
#'
#' @examples
#' \dontrun{
#' # Build taxonomy strings for DADA2 export
#' tax_strings <- rb_build_taxonomy_string(refs, style = "plain")
#'
#' # Standardize species names
#' valid_names <- rb_standardize_species(
#'   raw_names,
#'   reference_table = name_mapping
#' )
#'
#' # Join taxonomy onto curated data
#' final_data <- rb_join_taxonomy(
#'   qc_data,
#'   taxonomy_table = tax_table,
#'   defaults = list(kingdom = "Animalia", phylum = "Chordata")
#' )
#' }
#'
#' @name rb_taxonomy_utils
#' @family taxonomy
NULL

#' @rdname rb_taxonomy_utils
#' @export
rb_build_taxonomy_string <- function(data, style = c("plain", "rank_prefix")) {
  style <- match.arg(style)
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  rb_required_columns(data, ranks)
  prefixes <- c("k", "p", "c", "o", "f", "g", "s")
  apply(data[, ranks, drop = FALSE], 1, function(row) {
    vals <- as.character(row)
    vals[is.na(vals) | vals == ""] <- "unassigned"
    if (identical(style, "rank_prefix")) vals <- paste0(prefixes, "__", vals)
    paste(vals, collapse = ";")
  })
}

#' @rdname rb_taxonomy_utils
#' @export
rb_standardize_species <- function(names, reference_table,
                                    name_col = "alternative_clean",
                                    valid_col = "valid_name",
                                    taxonomy_lookup = NULL) {
  rb_required_columns(reference_table, c(name_col, valid_col))
  names_clean <- tolower(names)

  vapply(names_clean, function(nm) {
    match <- reference_table[[valid_col]][reference_table[[name_col]] == nm]
    match <- match[!is.na(match)]
    if (length(match) > 0) return(match[[1]])

    if (!is.null(taxonomy_lookup)) {
      looked_up <- tryCatch(taxonomy_lookup(nm), error = function(e) NA_character_)
      if (!is.na(looked_up)) return(looked_up)
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}

#' @rdname rb_taxonomy_utils
#' @export
rb_join_taxonomy <- function(data, taxonomy_table, by = "species",
                              defaults = list(), unassigned_label = "UNASSIGNED") {
  rb_required_columns(data, by)
  rb_required_columns(taxonomy_table, by)

  # Guard against merge()'s silent Cartesian join on duplicate keys.
  dup_species <- taxonomy_table[[by]][duplicated(taxonomy_table[[by]])]
  if (length(dup_species) > 0) {
    stop(
      "taxonomy_table has duplicate entries for: ", paste(unique(dup_species), collapse = ", "),
      ". merge() would silently multiply rows for these species -- ",
      "deduplicate taxonomy_table before calling rb_join_taxonomy().",
      call. = FALSE
    )
  }

  out <- merge(data, taxonomy_table, by = by, all.x = TRUE, suffixes = c("", ".tax"))

  rank_cols <- c("genus", "family", "order")
  for (col in rank_cols) {
    if (!col %in% names(out)) out[[col]] <- NA_character_
    out[[col]][is.na(out[[col]])] <- unassigned_label
  }

  for (rank_name in names(defaults)) {
    out[[rank_name]] <- defaults[[rank_name]]
  }

  out
}
