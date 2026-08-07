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
# DESTINATION: R/taxonomy.R (same file as rb_join_taxonomy from Module 5)
#
# Generalizes the species-matching-via-reference-table logic used
# identically across all 5 source blocks. check_cas() (a call out to
# Eschmeyer's Catalog of Fishes via rFishTaxa — fish-specific) becomes
# an optional `taxonomy_lookup` callback, NULL by default, so
# non-fish datasets skip that step entirely rather than needing a stub.
#
# NOTE: the original script's fuzzy-matching cross-check (comparing an
# NCBI header's extracted species name against the query species name
# via stringdist) is NOT included here — that's specific to how NCBI
# headers are parsed, not a general "standardize a name" concern. If
# you want that preserved, it belongs as an optional step inside
# rb_parse_source_table() for NCBI-style sources specifically, not in
# this function. Flag if you'd like me to add it there.
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

# --- Test ---
# ref <- data.frame(alternative_clean = "homo sapiens", valid_name = "Homo sapiens")
# rb_standardize_species("homo sapiens", reference_table = ref)
# expect: "Homo sapiens", no taxonomy_lookup call attempted
#
# rb_standardize_species("unknown species", reference_table = ref,
#                         taxonomy_lookup = function(nm) if (nm == "unknown species") "Resolved sp." else NA)
# expect: "Resolved sp." — confirms the fallback callback fires only
#         when the reference table lookup misses

 
#' MODULE 5> Join taxonomic lineage onto reference data with configurable defaults
#'
#' Generalizes Stage 3 of edna_ref_db_to_fasta.R. The original hardcoded
#' YZFishDB literals (kingdom = "Metazoa", phylum = "Chordata",
#' class = "Actinopteri", and the "Squalidus" genus-preservation special
#' case) are removed from the function body. Constant higher-rank values
#' now come from the caller via `defaults`; dataset-specific corrections
#' like the Squalidus case belong in the caller's own `taxonomy_table`
#' fix-up, not inside this function.
#'
#' @param data Data frame containing at least the `by` column.
#' @param taxonomy_table Reference table with columns: species, genus,
#'   family, order (kingdom/phylum/class are NOT expected here — supply
#'   those via `defaults`, since they're usually constant across an
#'   entire dataset, e.g. all sequences are Actinopterygii).
#' @param by Column to join on (default "species").
#' @param defaults Named list of constant higher-rank values applied to
#'   every row, e.g. list(kingdom = "Metazoa", phylum = "Chordata",
#'   class = "Actinopteri").
#' @param unassigned_label Label used when genus/family/order can't be
#'   resolved (default "UNASSIGNED", matching the original script).
rb_join_taxonomy <- function(data, taxonomy_table, by = "species",
                              defaults = list(), unassigned_label = "UNASSIGNED") {
  rb_required_columns(data, by)
  rb_required_columns(taxonomy_table, by)
 
  out <- merge(data, taxonomy_table, by = by, all.x = TRUE, suffixes = c("", ".tax"))
 
  rank_cols <- c("genus", "family", "order")
  for (col in rank_cols) {
    if (!col %in% names(out)) out[[col]] <- NA_character_
    out[[col]][is.na(out[[col]])] <- unassigned_label
  }
 
  # Caller-supplied constants replace the hardcoded
  # class = "Actinopteri" / phylum = "Chordata" / kingdom = "Metazoa"
  # from the original script.
  for (rank_name in names(defaults)) {
    out[[rank_name]] <- defaults[[rank_name]]
  }
 
  out
}
 
# --- Test ---
# tax <- data.frame(species = "Homo sapiens", genus = "Homo",
#                    family = "Hominidae", order = "Primates")
# rb_join_taxonomy(
#   data.frame(species = c("Homo sapiens", "Unknown sp.")),
#   taxonomy_table = tax,
#   defaults = list(kingdom = "Animalia", phylum = "Chordata", class = "Mammalia")
# )
# expect: row 1 fully resolved (genus/family/order from tax table);
#         row 2 genus/family/order = "UNASSIGNED";
#         BOTH rows get class = "Mammalia", not "Actinopteri"



