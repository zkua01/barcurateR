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
 
