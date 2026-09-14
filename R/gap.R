#' Barcode gap analysis for curated reference databases
#'
#' @description
#' Functions to assess whether reference sequences for a species form a
#' distinguishable barcode gap relative to sequences from other species.
#'
#' * `rb_run_barcode_gap()`: Runs barcode gap analysis across one or more
#'   markers. Aligns sequences with `DECIPHER::AlignSeqs()`, then computes
#'   per-species distance metrics.
#' * `rb_barcode_gap_species()`: Computes barcode gap metrics for a single
#'   species within an aligned marker dataset.
#' * `rb_diagnostic_sites()`: Counts fixed diagnostic nucleotide positions
#'   that distinguish focal sequences from non-focal sequences.
#'
#' @details
#' `rb_run_barcode_gap()` expects a standardized reference data frame with
#' at least species, marker/gene, and sequence columns. If a `qc_flag`
#' column is present, only sequences with `qc_flag == "pass"` are used.
#'
#' For each marker, the function:
#'
#' 1. Filters to the requested marker.
#' 2. Optionally removes alignment gap characters (`"-"`).
#' 3. Filters species with at least `min_seqs` sequences.
#' 4. Aligns sequences using `DECIPHER::AlignSeqs()`.
#' 5. Computes barcode gap metrics for each species.
#'
#' The returned data frame includes intra- and inter-specific distance
#' summaries, gap metrics, congeneric comparisons where available,
#' diagnostic site counts, sampling adequacy, taxonomic resolution,
#' assignment risk, and a recommended similarity threshold.
#'
#' If `parallel = TRUE`, a parallel cluster is used for the per-species
#' calculations. If cluster setup or execution fails, the function
#' automatically falls back to serial execution for that marker.
#'
#' Distance calculations use `rb_string_dist()`, which prefers
#' `pwalign::stringDist()` and falls back to `Biostrings::stringDist()`
#' on older Biostrings versions.
#'
#' @param data A standardized reference data frame.
#' @param markers Character vector of marker names to analyze.
#' @param species_col Name of the species column. Default `"species"`.
#' @param genus_col Name of the genus column. If `NULL`, genus is derived
#'   from the first token of the species name.
#' @param gene_col Name of the marker/gene column. Default `"seq_type"`.
#' @param sequence_col Name of the sequence column. Default `"sequence"`.
#' @param qc_flag_col Name of the QC flag column. Default `"qc_flag"`.
#'   If present, only sequences with `qc_flag == "pass"` are analyzed.
#' @param min_seqs Minimum number of sequences required per species.
#' @param parallel Logical. If `TRUE`, uses a parallel cluster where possible.
#' @param n_workers Optional number of parallel workers. If `NULL`, defaults
#'   to `parallel::detectCores() - 2`.
#' @param remove_gaps Logical. If `TRUE`, removes `"-"` characters from
#'   sequences before analysis. Default `TRUE`.
#' @param species_name Species to analyze in `rb_barcode_gap_species()`.
#' @param aligned_seqs Named character vector of aligned sequences.
#' @param seq_len Aligned sequence length used to normalize distances.
#' @param marker Marker name associated with the aligned sequences.
#' @param max_other Maximum number of non-focal sequences to use when
#'   computing global inter-specific distances. Default `2000`.
#' @param max_congeners Maximum number of congeneric sequences to use when
#'   computing congeneric distances. Default `500`.
#' @param focal_seqs Character vector of focal sequences for
#'   `rb_diagnostic_sites()`.
#' @param other_seqs Character vector of non-focal sequences for
#'   `rb_diagnostic_sites()`.
#' @param min_freq Minimum focal base frequency required for a position to
#'   be considered diagnostic. Default `0.9`.
#'
#' @return
#' * `rb_run_barcode_gap()` returns a data frame of barcode gap metrics,
#'   with one row per species per marker. Returns an empty data frame if
#'   no suitable data are available.
#' * `rb_barcode_gap_species()` returns a one-row data frame of barcode
#'   gap metrics, or `NULL` if the species does not meet the minimum
#'   sequence threshold or has no non-focal comparison sequences.
#' * `rb_diagnostic_sites()` returns an integer count of diagnostic sites.
#'
#' @section Output columns:
#' The data frame returned by `rb_run_barcode_gap()` and
#' `rb_barcode_gap_species()` includes:
#'
#' * `species`, `marker`, `n_sequences`
#' * `max_intra`, `min_inter`
#' * `gap_exists`, `gap_width`
#' * `p_overlap`, `q_overlap`
#' * `median_intra`, `median_inter`
#' * `intra_95`, `inter_5`, `percentile_gap`, `robust_gap`
#' * `n_congeners_in_db`, `min_congeneric`, `congeneric_gap`
#' * `congeneric_p_overlap`, `congeneric_q_overlap`
#' * `n_diagnostic_sites`
#' * `sampling_adequacy`, `taxonomic_resolution`
#' * `gap_may_be_undersampling_artifact`
#' * `error_msg`
#' * `assignment_risk`, `recommended_threshold`
#'
#' @examples
#' \dontrun{
#' gap_results <- rb_run_barcode_gap(
#'   final_data,
#'   markers = c("COI", "12S"),
#'   min_seqs = 5,
#'   parallel = FALSE
#' )
#'
#' head(gap_results)
#' }
#'
#' @name rb_barcode_gap_analysis
#' @family barcode gap
NULL

#' @rdname rb_barcode_gap_analysis
#' @export
rb_diagnostic_sites <- function(focal_seqs, other_seqs, min_freq = 0.9, max_other = 500) {
  if (length(other_seqs) > max_other) other_seqs <- sample(other_seqs, max_other)

  cm_f <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(focal_seqs), as.prob = TRUE)[c("A", "C", "G", "T"), ]
  cm_o <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(other_seqs), as.prob = TRUE)[c("A", "C", "G", "T"), ]

  pos_diag <- (cm_f["A", ] >= min_freq & cm_o["A", ] == 0) |
    (cm_f["C", ] >= min_freq & cm_o["C", ] == 0) |
    (cm_f["G", ] >= min_freq & cm_o["G", ] == 0) |
    (cm_f["T", ] >= min_freq & cm_o["T", ] == 0)
  sum(pos_diag, na.rm = TRUE)
}

#' Calculate pairwise sequence distances
#'
#' Internal helper that wraps `pwalign::stringDist()` where available,
#' falling back to older `Biostrings::stringDist()` versions.
#'
#' @keywords internal
#' @noRd
rb_string_dist <- function(x, method = "hamming") {
  if (requireNamespace("pwalign", quietly = TRUE)) {
    return(pwalign::stringDist(x, method = method))
  }
  
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    biostrings_version <- utils::packageVersion("Biostrings")
    
    if (biostrings_version < "2.77.1") {
      return(Biostrings::stringDist(x, method = method))
    }
  }
  
  stop(
    "Pairwise distance calculation requires the Bioconductor package 'pwalign'. ",
    "Install it with BiocManager::install('pwalign').",
    call. = FALSE
  )
}

#' @rdname rb_barcode_gap_analysis
#' @export
rb_barcode_gap_species <- function(species_name, data, aligned_seqs, seq_len,
                                    species_col = "species", genus_col = "genus",
                                    marker = NA_character_,
                                    min_seqs = 2,
                                    max_other = 2000, max_congeners = 500) {
  focal_idx <- which(data[[species_col]] == species_name)
  n_focal <- length(focal_idx)
  if (n_focal < min_seqs) return(NULL)

  if (is.na(marker) || !nzchar(marker)) {
    warning(
      "rb_barcode_gap_species() called with an unknown/blank marker for species '",
      species_name, "' -- skipping (barcode gap analysis requires a known marker).",
      call. = FALSE
    )
    return(NULL)
  }

  focal_seqs <- aligned_seqs[focal_idx]
  focal_genus <- data[[genus_col]][focal_idx[1]]

  other_idx <- which(data[[species_col]] != species_name)
  if (length(other_idx) == 0) return(NULL)
  if (length(other_idx) > max_other) other_idx <- sample(other_idx, max_other)
  other_seqs <- aligned_seqs[other_idx]

  congener_idx <- which(data[[genus_col]] == focal_genus & data[[species_col]] != species_name)
  congener_seqs <- if (length(congener_idx) > 0) {
    if (length(congener_idx) > max_congeners) congener_idx <- sample(congener_idx, max_congeners)
    aligned_seqs[congener_idx]
  } else {
    NULL
  }

  tryCatch({
    focal_dna <- Biostrings::DNAStringSet(focal_seqs)
    combined_global <- Biostrings::DNAStringSet(c(focal_seqs, other_seqs))

    intra_dist <- as.matrix(rb_string_dist(focal_dna, method = "hamming")) / seq_len
    max_intra <- max(intra_dist[upper.tri(intra_dist)], na.rm = TRUE)

    all_dist_global <- as.matrix(rb_string_dist(combined_global, method = "hamming")) / seq_len
    focal_other_dist <- all_dist_global[1:n_focal, (n_focal + 1):ncol(all_dist_global)]
    min_inter <- min(focal_other_dist, na.rm = TRUE)

    gap_exists <- max_intra < min_inter
    gap_width <- min_inter - max_intra

    median_intra <- stats::median(intra_dist[upper.tri(intra_dist)], na.rm = TRUE)
    median_inter <- stats::median(focal_other_dist, na.rm = TRUE)
    intra_95 <- stats::quantile(intra_dist[upper.tri(intra_dist)], 0.95, na.rm = TRUE)
    inter_5 <- stats::quantile(focal_other_dist, 0.05, na.rm = TRUE)
    percentile_gap <- inter_5 - intra_95

    if (!is.null(congener_seqs)) {
      combined_cong <- Biostrings::DNAStringSet(c(focal_seqs, congener_seqs))
      cong_dist <- as.matrix(rb_string_dist(combined_cong, method = "hamming")) / seq_len
      focal_cong_dist <- cong_dist[1:n_focal, (n_focal + 1):ncol(cong_dist)]
      min_congeneric <- min(focal_cong_dist, na.rm = TRUE)
      congeneric_gap <- min_congeneric - max_intra
      cong_p <- mean(intra_dist[upper.tri(intra_dist)] > min_congeneric, na.rm = TRUE)
      cong_q <- mean(focal_cong_dist < max_intra, na.rm = TRUE)
    } else {
      min_congeneric <- congeneric_gap <- cong_p <- cong_q <- NA_real_
    }

    n_diag <- rb_diagnostic_sites(focal_seqs, other_seqs)
    n_congeners <- length(unique(data[[species_col]][congener_idx]))

    sampling_adequacy <- dplyr::case_when(
      n_focal >= 20 & n_congeners >= 3 ~ "high",
      n_focal >= 10 & n_congeners >= 1 ~ "moderate",
      n_focal < 10 | n_congeners == 0 ~ "low"
    )
    taxonomic_resolution <- dplyr::case_when(
      gap_exists ~ "species",
      !is.na(congeneric_gap) & congeneric_gap > 0 ~ "species (congener-supported)",
      TRUE ~ "genus_or_higher"
    )

    data.frame(
      species = species_name, marker = marker, n_sequences = n_focal,
      max_intra = max_intra, min_inter = min_inter,
      gap_exists = gap_exists, gap_width = gap_width,
      p_overlap = round(mean(intra_dist[upper.tri(intra_dist)] > min_inter, na.rm = TRUE), 4),
      q_overlap = round(mean(focal_other_dist < max_intra, na.rm = TRUE), 4),
      median_intra = median_intra, median_inter = median_inter,
      intra_95 = intra_95, inter_5 = inter_5, percentile_gap = percentile_gap,
      robust_gap = percentile_gap > 0,
      n_congeners_in_db = n_congeners,
      min_congeneric = min_congeneric, congeneric_gap = congeneric_gap,
      congeneric_p_overlap = round(cong_p, 4), congeneric_q_overlap = round(cong_q, 4),
      n_diagnostic_sites = n_diag,
      sampling_adequacy = sampling_adequacy,
      taxonomic_resolution = taxonomic_resolution,
      gap_may_be_undersampling_artifact = (gap_exists & sampling_adequacy == "low"),
      error_msg = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      species = species_name, marker = marker, n_sequences = n_focal,
      max_intra = NA_real_, min_inter = NA_real_, gap_exists = NA,
      gap_width = NA_real_, p_overlap = NA_real_, q_overlap = NA_real_,
      median_intra = NA_real_, median_inter = NA_real_,
      intra_95 = NA_real_, inter_5 = NA_real_, percentile_gap = NA_real_, robust_gap = NA,
      n_congeners_in_db = NA_integer_, min_congeneric = NA_real_,
      congeneric_gap = NA_real_, congeneric_p_overlap = NA_real_, congeneric_q_overlap = NA_real_,
      n_diagnostic_sites = NA_integer_, sampling_adequacy = NA_character_,
      taxonomic_resolution = NA_character_, gap_may_be_undersampling_artifact = NA,
      error_msg = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

#' @rdname rb_barcode_gap_analysis
#' @export
rb_run_barcode_gap <- function(data, markers = c("COI", "12S", "16S", "genome"),
                                species_col = "species", genus_col = NULL,
                                gene_col = "seq_type", sequence_col = "sequence",
                                qc_flag_col = "qc_flag", min_seqs = 5,
                                parallel = TRUE, n_workers = NULL, remove_gaps = TRUE) {
  rb_required_columns(data, c(species_col, gene_col, sequence_col))

  if (is.null(genus_col) || !genus_col %in% names(data)) {
    data$genus <- sapply(strsplit(data[[species_col]], "\\s+"), `[`, 1)
    genus_col <- "genus"
  }

  process_one_marker <- function(marker) {
    message("Processing marker: ", marker)

    db <- data[data[[gene_col]] == marker, , drop = FALSE]
    if (!is.null(qc_flag_col) && qc_flag_col %in% names(db)) {
      db <- db[db[[qc_flag_col]] == "pass", , drop = FALSE]
    }
    if (isTRUE(remove_gaps)){
      db[[sequence_col]] <- gsub("-", "", db[[sequence_col]], fixed = TRUE)
    }
    species_counts <- table(db[[species_col]])
    db <- db[db[[species_col]] %in% names(species_counts[species_counts >= min_seqs]), , drop = FALSE]

    if (nrow(db) == 0) {
      message("  No data with n >= ", min_seqs, " for ", marker)
      return(NULL)
    }

    species_list <- unique(db[[species_col]])
    seqs <- Biostrings::DNAStringSet(db[[sequence_col]])
    names(seqs) <- if ("unique_code" %in% names(db)) db$unique_code else paste0("seq_", seq_len(nrow(db)))

    aligned <- DECIPHER::AlignSeqs(seqs, verbose = FALSE)
    aligned_chars <- stats::setNames(as.character(aligned), names(seqs))
    seq_len_val <- Biostrings::width(aligned)[1]

    gap_list <- if (isTRUE(parallel)) {
      workers <- if (is.null(n_workers)) max(1, parallel::detectCores() - 2) else n_workers
      tryCatch({
        cl <- parallel::makeCluster(workers)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterExport(cl, varlist = "rb_diagnostic_sites")
        pbapply::pblapply(
          species_list, rb_barcode_gap_species,
          data = db, aligned_seqs = aligned_chars, seq_len = seq_len_val,
          species_col = species_col, genus_col = genus_col, marker = marker,
          min_seqs = min_seqs, cl = cl
        )
      }, error = function(e) {
        warning(
          "Parallel barcode gap analysis failed for marker '", marker, "' (",
          conditionMessage(e), ") -- falling back to serial execution.",
          call. = FALSE
        )
        lapply(
          species_list, rb_barcode_gap_species,
          data = db, aligned_seqs = aligned_chars, seq_len = seq_len_val,
          species_col = species_col, genus_col = genus_col, marker = marker,
          min_seqs = min_seqs
        )
      })
    } else {
      lapply(
        species_list, rb_barcode_gap_species,
        data = db, aligned_seqs = aligned_chars, seq_len = seq_len_val,
        species_col = species_col, genus_col = genus_col, marker = marker,
        min_seqs = min_seqs
      )
    }

    valid_list <- gap_list[!vapply(gap_list, is.null, logical(1))]
    if (length(valid_list) == 0) return(NULL)
    gap_res <- do.call(rbind, valid_list)

    success_mask <- is.na(gap_res$error_msg)
    thresh_95 <- if (any(success_mask)) {
      stats::quantile(gap_res$max_intra[success_mask], 0.95, na.rm = TRUE)
    } else {
      NA_real_
    }

    gap_res$assignment_risk <- dplyr::case_when(
      !is.na(gap_res$error_msg) ~ "failed",
      gap_res$gap_exists & gap_res$gap_width > 0.02 ~ "low",
      gap_res$gap_exists & gap_res$gap_width <= 0.02 ~ "moderate",
      TRUE ~ "high"
    )
    gap_res$recommended_threshold <- ifelse(success_mask, round(max(thresh_95 * 1.05, 0.97), 3), NA_real_)
    gap_res
  }

  results <- lapply(markers, process_one_marker)
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) return(data.frame())
  do.call(rbind, results)
}
