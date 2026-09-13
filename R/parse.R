# ============================================================
# R/parse.R
# Source parsing: marker extraction, per-source table parsing,
# primer trimming, combining sources, ambiguous-sequence resolution.
# ============================================================


# ------------------------------------------------------------
# rb_default_marker_patterns() / rb_extract_markers()
#
# Generalizes extract_markers(). marker_patterns, genome_pattern,
# multi_gene_pattern, and other_gene_pattern are now parameters
# instead of literals baked into the function body — a dataset with a
# different marker set (e.g. ITS/rbcL for plants, or a 5-marker fish
# panel) supplies its own marker_patterns list instead of editing the
# function.
# ------------------------------------------------------------

rb_default_marker_patterns <- function() {
  list(
    "12S" = "12s|mt-rnr1|12s ribosomal rna",
    "16S" = "16s|mt-rnr2|16s ribosomal rna",
    "COI" = "coi|co1|cox1|cox i|cox1|cytochrome c oxidase subunit i|cytochrome c oxidase subunit i|cytochrome oxidase subunit i",
    "CYTB" = "cytb|cytochrome b|mt-cytb",
    "18S" = "18s|18s ribosomal rna|18s rrna"
  )
}

#' rb_extract_markers
#' Extract marker type + completeness from a sequence description
#' @param description Free-text sequence description (e.g. a GenBank
#'   definition line).
#' @param marker_patterns Named list mapping marker name -> detection
#'   regex. Default covers 12S/16S/COI; supply your own for other
#'   marker sets (e.g. list(ITS = "its")).
#' @param genome_pattern Regex identifying a complete-genome record.
#' @param multi_gene_pattern Regex identifying a multi-gene record.
#' @param other_gene_pattern Regex identifying non-marker gene mentions
#'   (tRNA, cytb, ND, ATP genes) used to flag an "other" component in a
#'   multi-gene record.

rb_extract_markers <- function(description,
                                marker_patterns = rb_default_marker_patterns(),
                                genome_pattern = "complete genome|mitochondrion, complete genome|mitogenome",
                                multi_gene_pattern = paste0(
                                  "genes for|gene for|and|;|,.*and|partial.*complete|",
                                  "complete.*partial|complete.*complete|rrna.*rrna|",
                                  "rrna.*trna|trna.*rrna"
                                ),
                                other_gene_pattern = "trna|transfer rna|nd[0-9]|atp[0-9]|nad[0-9]") {
  if(length(description) != 1 || is.na(description) || !nzchar(description)){
    return("other_unknown")
  }
  
  desc_lower <- tolower(description)

  if (stringr::str_detect(desc_lower, genome_pattern)) {
    return("genome_complete")
  }

  detect_completeness <- function(marker_pattern) {
    p <- paste0("(?:", marker_pattern, ")")
    
    partial_regex <- paste0(p, ".*partial|partial.*", p)
    complete_regex <- paste0(p, ".*complete|complete.*", p)
    
    if (stringr::str_detect(desc_lower, partial_regex)) {
      "partial"
    } else if (stringr::str_detect(desc_lower, complete_regex)) {
      "complete"
    } else {
      "unknown"
    }
  }

  is_multi_gene <- stringr::str_detect(desc_lower, multi_gene_pattern)
  
  present_markers <- names(marker_patterns)[
    vapply(marker_patterns, function(p) stringr::str_detect(desc_lower, p), logical(1))
  ]

    if (is_multi_gene) {
    markers <- character(0)
    for (marker_name in present_markers) {
      comp <- detect_completeness(marker_patterns[[marker_name]])
      markers <- c(markers, paste0(marker_name, "_", comp))
    }
    has_other <- stringr::str_detect(desc_lower, other_gene_pattern)
    if (has_other) {
      comp_other <- if (stringr::str_detect(desc_lower, "partial")) "partial"
                    else if (stringr::str_detect(desc_lower, "complete")) "complete"
                    else "unknown"
      markers <- c(markers, paste0("other_", comp_other))
    }
    if (length(markers) > 0) return(paste(unique(markers), collapse = ";"))
    return("other_unknown")
  }

  if (length(present_markers) > 0) {
    marker_name <- present_markers[[1]]
    comp <- detect_completeness(marker_patterns[[marker_name]])
    return(paste0(marker_name, "_", comp))
  }

  "other_unknown"
}


# ------------------------------------------------------------
# rb_parse_source_table()
#' Parse one source's raw sequence table into the standard schema
#'
#' @param raw Raw data.frame as read from a source file.
#' @param source_name Short label for this source (e.g. "bold", "ncbi").
#' @param column_map Named character vector: standard name -> column
#'   name in `raw`. Passed to rb_standardize_columns() (R/utils.R).
#' @param marker_fn Function to derive seq_type from a description
#'   column, default rb_extract_markers().
#' @param description_col Optional name (in `raw`, before renaming) of
#'   a free-text description column to run marker_fn over. If omitted
#'   and no seq_type column exists after renaming, seq_type defaults
#'   to "other_unknown".
#
# Generalizes the 5 separate hardcoded per-source blocks (NCBI/BOLD/
# MitoFish/MIDORI2/local) into one function driven by a column_map,
# using rb_standardize_columns() from Module 5 to do the renaming.
# ------------------------------------------------------------

rb_parse_source_table <- function(raw, source_name, column_map,
                                   marker_fn = rb_extract_markers,
                                   description_col = NULL) {
  parsed <- rb_standardize_columns(raw, column_map)
  parsed$source <- source_name
  parsed$sequence <- rb_clean_sequence(parsed$sequence)
  parsed$length <- nchar(parsed$sequence)

  if (!is.null(description_col) && description_col %in% names(raw)) {
    parsed$seq_type <- vapply(raw[[description_col]], marker_fn, character(1))
  } else if (!"seq_type" %in% names(parsed)) {
    parsed$seq_type <- "other_unknown"
  }

  parsed
}

# ------------------------------------------------------------
# rb_reverse_complement() / rb_remove_primers()
# Already dataset-agnostic logic — pulled out of the BOLD-only code
# path so any source with primer-flanked sequences can call it.
# ------------------------------------------------------------

rb_reverse_complement <- function(dna_seq) {
  complement <- chartr("ATCG", "TAGC", dna_seq)
  stringi::stri_reverse(complement)
}
           
#' Trim flanking primers off a sequence, if present
rb_remove_primers <- function(sequence, f_primer_seqs = NA_character_,
                               r_primer_seqs = NA_character_) {
  if (all(is.na(f_primer_seqs)) && all(is.na(r_primer_seqs))) return(sequence)

  for (f_primer in f_primer_seqs) {
    if (!is.na(f_primer) && startsWith(sequence, f_primer)) {
      sequence <- substr(sequence, nchar(f_primer) + 1, nchar(sequence))
      break
    }
  }
  for (r_primer in r_primer_seqs) {
    if (!is.na(r_primer)) {
      rc_primer <- rb_reverse_complement(r_primer)
      if (endsWith(sequence, rc_primer)) {
        sequence <- substr(sequence, 1, nchar(sequence) - nchar(rc_primer))
        break
      }
    }
  }
  sequence
}


# ------------------------------------------------------------
# rb_combine_sources()
 #' Combine multiple standardized source tables, deduplicating by
#' species + sequence.
# ------------------------------------------------------------

rb_combine_sources <- function(source_list, species_col = "species", sequence_col = "sequence") {
  combined <- dplyr::bind_rows(source_list)
  combined[[sequence_col]] <- rb_clean_sequence(combined[[sequence_col]])
  combined <- combined[!duplicated(combined[c(species_col, sequence_col)]), , drop = FALSE]
  row.names(combined) <- NULL
  combined
}


# ------------------------------------------------------------
# rb_resolve_ambiguous()
#' Resolve sequences that were assigned to multiple species
#'
#' @param resolution_table Table with an `action` column ("keep" or
#'   "combine") resolving each ambiguous sequence.
#' @param combine_cols Columns that may legitimately differ across an
#'   ambiguous group and should be pipe-joined when combined (e.g.
#'   species/source/sequence_id/seq_type). Every other column is taken
#'   via first() on the assumption it's identical across the group
#'   (true by construction for e.g. `length`, since the group shares
#'   one sequence).
#' @param on_unresolved What to do when ambiguous sequences exist and
#'   no resolution_table is supplied: "stop" (default, matches
#'   original script), "warn" (return data unchanged), or "drop"
#'   (remove the ambiguous rows).
#' @param flag_output_path Optional path to write the ambiguous rows
#'   for manual review before resolving (only used in the
#'   no-resolution-table case).
# ------------------------------------------------------------

rb_resolve_ambiguous <- function(data, species_col = "species", sequence_col = "sequence",
                                  resolution_table = NULL,
                                  on_unresolved = "stop",
                                  flag_output_path = NULL) {
  on_unresolved <- rb_match_arg(on_unresolved, c("stop","warn","drop"))

  dup_check <- data %>%
    dplyr::group_by(.data[[sequence_col]]) %>%
    dplyr::filter(dplyr::n_distinct(.data[[species_col]]) > 1) %>%
    dplyr::ungroup()

  if (nrow(dup_check) == 0) {
    return(data)
  }

  if (is.null(resolution_table)) {
    if (!is.null(flag_output_path)) {
      utils::write.csv(dup_check, flag_output_path, row.names = FALSE)
    }
    msg <- sprintf(
      "%d ambiguous sequences found (assigned to multiple species) with no resolution_table supplied.",
      nrow(dup_check)
    )
    if (on_unresolved == "stop") stop(msg, call. = FALSE)
    if (on_unresolved == "warn") { warning(msg, call. = FALSE); return(data) }
    if (on_unresolved == "drop") {
      warning(paste(msg, "Dropping them."), call. = FALSE)
      return(dplyr::anti_join(data, dup_check, by = sequence_col))
    }
  }

  # resolution_table expected to carry an `action` column: "keep" or "combine"
  if(!"action" %in% names(resolution_table)){
    stop("resolution_table must contain an 'action' column.", call.=FALSE)
  }
  action_normalized <- tolower(trimws(as.character(resolution_table$action)))
  
  action_normalized[is.na(action_normalized) | action_normalized == ""] <- "drop"
  action_normalized[!action_normalized %in% c("keep","combine","drop")] <- "drop"
  resolution_table$action <- action_normalized
  
  action_counts <- table(
    factor(resolution_table$action, levels = c("keep", "combine", "drop"))
  )
  
  if (action_counts[["drop"]] > 0) {
    warning(
      sprintf(
        "Dropped %d resolution row(s) because action was missing, blank, or 'drop'. Action counts: keep = %d, combine = %d, drop/blank = %d.",
        action_counts[["drop"]],
        action_counts[["keep"]],
        action_counts[["combine"]],
        action_counts[["drop"]]
      ),
      call. = FALSE
    )
  } else {
    message(
      sprintf(
        "Resolution table action counts: keep = %d, combine = %d.",
        action_counts[["keep"]],
        action_counts[["combine"]]
      )
    )
  }
  
  keep_rows <- resolution_table[resolution_table$action == "keep", , drop = FALSE]
  keep_rows$action <- NULL

  combine_rows <- resolution_table[resolution_table$action == "combine", , drop = FALSE]
  if (nrow(combine_rows) > 0) {
    combine_rows <- combine_rows %>%
      dplyr::group_by(.data[[sequence_col]]) %>%
      dplyr::summarise(
        dplyr::across(-dplyr::all_of("action"), ~ paste(unique(.x), collapse = "|")),
        .groups = "drop"
      )
  }

  data %>%
    dplyr::anti_join(resolution_table, by = sequence_col) %>%
    dplyr::bind_rows(keep_rows) %>%
    dplyr::bind_rows(combine_rows) %>%
    dplyr::distinct(.data[[species_col]], .data[[sequence_col]], .keep_all = TRUE)
}



