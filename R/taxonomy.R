# ============================================================
# R/taxonomy.R
# Species-name standardization and taxonomic lineage joining.
# (rb_build_taxonomy_string(), if it already lives in this file, is
# unchanged by this refactor and should stay here alongside these.)
# ============================================================

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


# ------------------------------------------------------------
# rb_standardize_species()
#
#' Standardize a species name against a reference table, with an
#' optional external taxonomy lookup as fallback
#'
#' @param reference_table Table mapping known alternative names to the
#'   currently valid name.
#' @param name_col Column in reference_table holding alternative/known
#'   names (lowercased for matching).
#' @param valid_col Column in reference_table holding the valid name.
#' @param taxonomy_lookup Optional callback `function(name) -> valid_name`
#'   used when reference_table has no match (e.g. an external catalog
#'   API). NULL by default, so non-fish datasets skip this step
#'   entirely rather than needing a stub implementation.
#'
#' Note: the original script's fuzzy-match cross-check (comparing an
#' NCBI header's extracted species name against the query species name
#' via stringdist) is intentionally NOT included here — it's specific
#' to NCBI header parsing, not general name standardization. If wanted,
#' it belongs inside rb_parse_source_table() (R/parse.R) for
#' NCBI-style sources.
# ------------------------------------------------------------

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

 
#' Join taxonomic lineage onto reference data with configurable defaults
#'
#' @param taxonomy_table Reference table with columns: species (or
#'   whatever `by` is), genus, family, order.
#' @param by Column to join on (default "species").
#' @param defaults Named list of constant higher-rank values applied
#'   to every row, e.g. list(kingdom = "Metazoa", phylum = "Chordata",
#'   class = "Actinopteri").
#' @param unassigned_label Label used when genus/family/order can't be
#'   resolved (default "UNASSIGNED").
rb_join_taxonomy <- function(data, taxonomy_table, by = "species",
                              defaults = list(), unassigned_label = "UNASSIGNED") {
  rb_required_columns(data, by)
  rb_required_columns(taxonomy_table, by)

  # Guard against merge()'s silent Cartesian join on duplicate keys.
  dup_species <- taxonomy_table[[by]][duplicated(taxonomy_table[[by]])]
  if (length(dup_species) > 0) {
    stop(
      "taxonomy_table has duplicate entries for: ", paste(unique(dup_species), collapse = ", "),
      ". merge() would silently multiply rows for these species — ",
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


