#' Parse and standardize source sequence tables
#'
#' @description
#' Functions to extract marker information from sequence descriptions and
#' parse raw source tables into the standardized barcurateR schema.
#'
#' * `rb_default_marker_patterns()`: Returns the default regular expressions
#'   used to detect common mitochondrial and ribosomal markers.
#' * `rb_extract_markers()`: Extracts marker type and completeness from a
#'   free-text sequence description.
#' * `rb_parse_source_table()`: Parses a raw source table into the
#'   standardized schema used by downstream QC and curation functions.
#'
#' @details
#' The standardized schema used by barcurateR expects at least:
#'
#' * `species`
#' * `sequence`
#' * `seq_type`
#' * `source`
#'
#' Additional columns such as `sequence_id`, `description`, `length`,
#' and source-specific metadata may also be retained.
#'
#' `rb_parse_source_table()` uses a user-supplied `column_map` to rename
#' raw columns into the standardized schema. If a description column is
#' provided, marker type is inferred using `rb_extract_markers()` and then
#' simplified with the internal helper `rb_simplify_seq_type()`.
#'
#' @param description Free-text sequence description, for example a GenBank
#'   definition line.
#' @param marker_patterns Named list mapping marker names to detection
#'   regular expressions.
#' @param genome_pattern Regular expression identifying complete genome
#'   records.
#' @param multi_gene_pattern Regular expression identifying multi-gene
#'   records.
#' @param other_gene_pattern Regular expression identifying non-marker gene
#'   mentions such as tRNA, ND, or ATP genes.
#' @param raw Raw data frame as read from a source file.
#' @param source_name Short label for the source, for example `"ncbi"` or
#'   `"bold"`.
#' @param column_map Named character vector mapping standardized column
#'   names to column names in `raw`.
#' @param marker_fn Function used to derive marker type from a description
#'   column. Default is `rb_extract_markers()`.
#' @param description_col Optional name of a free-text description column
#'   in `raw`.
#' @param min_length Minimum sequence length to retain. Default `1`.
#'
#' @return
#' * `rb_default_marker_patterns()` returns a named list of regular
#'   expressions.
#' * `rb_extract_markers()` returns a character string describing the
#'   detected marker and completeness, for example `"COI_partial"`,
#'   `"12S_unknown"`, `"genome_complete"`, `"multi_marker"`, or
#'   `"other_unknown"`.
#' * `rb_parse_source_table()` returns a standardized data frame.
#'
#' @examples
#' \dontrun{
#' raw_ncbi <- read.csv("ncbi_sequences.csv", stringsAsFactors = FALSE)
#'
#' parsed_ncbi <- rb_parse_source_table(
#'   raw_ncbi,
#'   source_name = "ncbi",
#'   column_map = c(
#'     sequence_id = "sequence_id",
#'     species = "species_query",
#'     sequence = "sequence"
#'   ),
#'   description_col = "description"
#' )
#' }
#'
#' @name rb_parse_sources
#' @family parsing
NULL

#' @rdname rb_parse_sources
#' @export
rb_default_marker_patterns <- function() {
  list(
    "12S" = "12s|mt-rnr1|12s ribosomal rna",
    "16S" = "16s|mt-rnr2|16s ribosomal rna",
    "COI" = "coi|co1|cox1|cox i|cox1|cytochrome c oxidase subunit i|cytochrome c oxidase subunit i|cytochrome oxidase subunit i",
    "CYTB" = "cytb|cytochrome b|mt-cytb",
    "18S" = "18s|18s ribosomal rna|18s rrna"
  )
}

#' Simplify marker labels for downstream QC
#'
#' Internal helper that converts detailed marker labels returned by
#' `rb_extract_markers()` into clean `seq_type` values used by QC checks.
#'
#' @keywords internal
#' @noRd
rb_simplify_seq_type <- function(marker_label) {
  vapply(marker_label, function(label) {
    
    # Handle empty / NA
    if (is.na(label) || !nzchar(label)) {
      return("other_unknown")
    }
    
    # Handle multi-gene records (contain semicolons)
    if (grepl(";", label)) {
      return("multi_marker")
    }
    
    # Standardize genome_complete
    if (grepl("genome_complete|genome", label, ignore.case = TRUE)) {
      return("genome_complete")
    }
    
    # Strip completeness suffix (_partial, _complete, _unknown)
    marker_base <- sub("_(partial|complete|unknown)$", "", label)
    
    # Handle other_* patterns
    if (grepl("^other", marker_base, ignore.case = TRUE)) {
      return("other_unknown")
    }
    
    # Uppercase for matching
    marker_upper <- toupper(marker_base)
    
    # Check against known markers
    known_markers <- c("COI", "12S", "16S", "18S", "CYTB")
    
    if (marker_upper %in% known_markers) {
      return(marker_upper)
    }
    
    "other_unknown"
  }, character(1), USE.NAMES = FALSE)
}

#' @rdname rb_parse_sources
#' @export
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

#' @rdname rb_parse_sources
#' @export
rb_parse_source_table <- function(raw, source_name, column_map,
                                  marker_fn = rb_extract_markers,
                                  description_col = NULL,
                                  min_length = 1) {
  parsed <- rb_standardize_columns(raw, column_map)
  parsed$source <- source_name
  parsed$sequence <- rb_clean_sequence(parsed$sequence)
  parsed$length <- nchar(parsed$sequence)
  
  if (!is.null(description_col) && description_col %in% names(raw)) {
    parsed$seq_type <- vapply(raw[[description_col]], marker_fn, character(1))
    parsed$seq_type <- rb_simplify_seq_type(parsed$seq_type)
  } else if (!"seq_type" %in% names(parsed)) {
    parsed$seq_type <- "other_unknown"
  }
  
  # Remove empty or zero-length sequences
  parsed <- parsed[parsed$length >= min_length, , drop = FALSE]
  row.names(parsed) <- NULL
  
  parsed
}

#' Prepare and trim raw DNA sequences
#'
#' @description
#' Utility functions for preparing raw DNA sequences during source
#' standardization.
#'
#' * `rb_reverse_complement()`: Returns the reverse complement of a DNA
#'   sequence.
#' * `rb_remove_primers()`: Removes matching forward and reverse primer
#'   sequences from the ends of a sequence, if present.
#'
#' @param dna_seq Character vector of DNA sequences.
#' @param sequence Character vector of DNA sequences to trim.
#' @param f_primer_seqs Character vector of forward primer sequences.
#'   Use `NA` if no forward primers should be removed.
#' @param r_primer_seqs Character vector of reverse primer sequences.
#'   Use `NA` if no reverse primers should be removed.
#'
#' @return
#' A character vector of transformed sequences.
#'
#' @details
#' `rb_remove_primers()` performs simple exact matching at sequence ends.
#' Forward primers are removed from the start of the sequence. Reverse
#' primers are reverse-complemented and removed from the end of the
#' sequence.
#'
#' @examples
#' \dontrun{
#' rb_reverse_complement("ATGC")
#'
#' rb_remove_primers(
#'   "ACGTACGTACGT",
#'   f_primer_seqs = "ACGT",
#'   r_primer_seqs = NA
#' )
#' }
#'
#' @name rb_sequence_prep
#' @family parsing
NULL

#' @rdname rb_sequence_prep
#' @export
rb_reverse_complement <- function(dna_seq) {
  complement <- chartr("ATCG", "TAGC", dna_seq)
  stringi::stri_reverse(complement)
}
           
#' @rdname rb_sequence_prep
#' @export
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

#' Combine parsed sources and resolve ambiguous sequences
#'
#' @description
#' Functions to combine multiple standardized source tables and resolve
#' sequences that are assigned to multiple species.
#'
#' * `rb_combine_sources()`: Combines parsed source tables and removes
#'   duplicated species-sequence combinations.
#' * `rb_resolve_ambiguous()`: Detects and resolves sequences assigned to
#'   multiple species.
#'
#' @param source_list List of standardized source data frames.
#' @param species_col Name of the species column. Default `"species"`.
#' @param sequence_col Name of the sequence column. Default `"sequence"`.
#' @param data Standardized data frame to check for ambiguous sequences.
#' @param resolution_table Optional data frame resolving ambiguous
#'   sequences. Must contain an `action` column with values `"keep"`,
#'   `"combine"`, or `"drop"`.
#' @param on_unresolved Action to take when ambiguous sequences are found
#'   and no `resolution_table` is supplied. One of `"stop"`, `"warn"`,
#'   or `"drop"`. Default `"stop"`.
#' @param flag_output_path Optional path to write ambiguous rows for manual
#'   review before resolving.
#'
#' @return
#' * `rb_combine_sources()` returns a combined data frame with duplicated
#'   species-sequence combinations removed.
#' * `rb_resolve_ambiguous()` returns a resolved data frame.
#'
#' @details
#' Ambiguous sequences are sequences that appear under more than one
#' species. These must be resolved before QC because downstream QC and
#' curation assume that each sequence is assigned to a single species.
#'
#' If a `resolution_table` is supplied, it should contain one row per
#' ambiguous sequence and an `action` column:
#'
#' * `"keep"`: keep the specified resolution.
#' * `"combine"`: combine records across species, pipe-joining values that
#'   differ across the ambiguous group.
#' * `"drop"`: remove the ambiguous sequence.
#'
#' If no `resolution_table` is supplied, `on_unresolved` controls the
#' behavior:
#'
#' * `"stop"`: stop with an error.
#' * `"warn"`: warn and return the data unchanged.
#' * `"drop"`: warn and remove ambiguous rows.
#'
#' @examples
#' \dontrun{
#' combined <- rb_combine_sources(
#'   list(ncbi_parsed, bold_parsed),
#'   species_col = "species",
#'   sequence_col = "sequence"
#' )
#'
#' resolved <- rb_resolve_ambiguous(
#'   combined,
#'   on_unresolved = "drop"
#' )
#' }
#'
#' @name rb_combine_resolve
#' @family parsing
NULL

#' @rdname rb_combine_resolve
#' @export
rb_combine_sources <- function(source_list, species_col = "species", sequence_col = "sequence") {
  combined <- dplyr::bind_rows(source_list)
  combined[[sequence_col]] <- rb_clean_sequence(combined[[sequence_col]])
  combined <- combined[!duplicated(combined[c(species_col, sequence_col)]), , drop = FALSE]
  row.names(combined) <- NULL
  combined
}

#' @rdname rb_combine_resolve
#' @export
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
