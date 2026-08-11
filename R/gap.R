
# ============================================================
# MODULE 6 — Barcode gap analysis (replaces the standalone gap script)
# DESTINATION: R/gap.R (new file)
#
# Output column names here are deliberately kept identical to the
# original script's (species, marker, gap_exists, gap_width, etc.) so
# results feed directly into your EXISTING rb_barcode_gap() /
# rb_read_barcode_gap() reporting functions without any schema
# translation step.
# ============================================================
 
 
# ------------------------------------------------------------
# rb_diagnostic_sites()
# Already dataset-agnostic — extracted as-is, with the previously
# hardcoded `500`-sequence subsampling cap now a parameter.
# ------------------------------------------------------------
 
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
 
# --- Test ---
# focal <- c("AAAA","AAAA","AAAA")
# other <- c("TTTT","TTTT","TTTT")
# rb_diagnostic_sites(focal, other, min_freq = 0.9)
# expect: 4 (every position is diagnostic)
 
 
# ------------------------------------------------------------
# rb_barcode_gap_species()
# Generalizes compute_gap(). `species_col`/`genus_col` were previously
# hardcoded as `ref_df$species`/`ref_df$genus`; `max_other` (2000) and
# `max_congeners` (500) subsampling caps are now parameters; `marker`
# is passed explicitly instead of being captured from an enclosing
# scope variable.
# ------------------------------------------------------------
 
rb_barcode_gap_species <- function(species_name, data, aligned_seqs, seq_len,
                                    species_col = "species", genus_col = "genus",
                                    marker = NA_character_,
                                    max_other = 2000, max_congeners = 500) {
  focal_idx <- which(data[[species_col]] == species_name)
  n_focal <- length(focal_idx)
  if (n_focal < 2) return(NULL)
 
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
 
    intra_dist <- as.matrix(Biostrings::stringDist(focal_dna, method = "hamming")) / seq_len
    max_intra <- max(intra_dist[upper.tri(intra_dist)], na.rm = TRUE)
 
    all_dist_global <- as.matrix(Biostrings::stringDist(combined_global, method = "hamming")) / seq_len
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
      cong_dist <- as.matrix(Biostrings::stringDist(combined_cong, method = "hamming")) / seq_len
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
 
# --- Test ---
# toy <- data.frame(species = c("sp1","sp1","sp1","sp2","sp2","sp2"),
#                    genus = c("Genus1","Genus1","Genus1","Genus2","Genus2","Genus2"))
# aligned <- stats::setNames(
#   c("AAAAAAAA","AAAAAAAT","AAAAAAAA","TTTTTTTT","TTTTTTTA","TTTTTTTT"),
#   paste0("seq_", 1:6)
# )
# rb_barcode_gap_species("sp1", toy, aligned_seqs = aligned, seq_len = 8,
#                         genus_col = "genus", marker = "COI")
# expect: a 1-row data.frame with gap_exists = TRUE (sp1 and sp2 are
#         clearly distinguishable), n_sequences = 3, error_msg = NA
 
 
# ------------------------------------------------------------
# rb_run_barcode_gap()
# Generalizes process_marker() + the top-level purrr::map(markers, ...)
# driver. `markers` is now a parameter instead of a hardcoded
# c("COI","12S","16S","genome") vector, so a 2-marker-only study (say,
# COI + 12S) doesn't need to touch the function body. Column names for
# species/genus/gene/sequence/qc_flag are all parameterized rather than
# assumed. If no genus column is supplied, genus is derived from the
# first token of the species name (matching the original script's
# `sapply(strsplit(species, "\\s+"), \`[\`, 1)` fallback).
# ------------------------------------------------------------
 
rb_run_barcode_gap <- function(data, markers = c("COI", "12S", "16S", "genome"),
                                species_col = "species", genus_col = NULL,
                                gene_col = "seq_type", sequence_col = "sequence",
                                qc_flag_col = "qc_flag", min_seqs = 5,
                                parallel = TRUE, n_workers = NULL) {
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
 
    if (isTRUE(parallel)) {
      workers <- if (is.null(n_workers)) max(1, parallel::detectCores() - 2) else n_workers
      cl <- parallel::makeCluster(workers)
      parallel::clusterExport(cl, varlist = "rb_diagnostic_sites")
      on.exit(parallel::stopCluster(cl), add = TRUE)
      gap_list <- pbapply::pblapply(
        species_list, rb_barcode_gap_species,
        data = db, aligned_seqs = aligned_chars, seq_len = seq_len_val,
        species_col = species_col, genus_col = genus_col, marker = marker,
        cl = cl
      )
    } else {
      gap_list <- lapply(
        species_list, rb_barcode_gap_species,
        data = db, aligned_seqs = aligned_chars, seq_len = seq_len_val,
        species_col = species_col, genus_col = genus_col, marker = marker
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
 
# --- Test ---
# toy <- data.frame(
#   species = rep(c("sp1", "sp2"), each = 6),
#   seq_type = "COI", qc_flag = "pass",
#   sequence = c(rep("ACGTACGTACGT", 6), rep("TTTTGGGGCCCC", 6))
# )
# result <- rb_run_barcode_gap(toy, markers = "COI", min_seqs = 5, parallel = FALSE)
# expect: runs without requiring 12S/16S/genome columns to exist at all,
#         since markers = "COI" only; returns 2 rows (sp1, sp2) with
#         gap_exists = TRUE for both
#
# Also confirms genus derivation works without an explicit genus column:
# unique(sapply(strsplit(toy$species, "\\s+"), `[`, 1))
# expect: c("sp1", "sp2") used as genus fallback (single-token species
#         names here just mean genus == species, which is fine for the test)
