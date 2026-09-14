#' Curate a DNA reference sequence database
#'
#' @description
#' `rb_curate_reference()` is the main user-facing wrapper that orchestrates
#' the full reference database curation pipeline:
#'
#' 1. First-pass quality control (contaminants, NUMTs, codons, rRNA).
#' 2. Optional machine-learning classification of unknown sequences.
#' 3. Second-pass quality control using updated marker labels.
#' 4. Optional taxonomic lineage joining and validation.
#' 5. Optional barcode gap analysis.
#' 6. Optional SQLite database export.
#'
#' @param data A standardized data frame containing at least `species`,
#'   `sequence`, and `seq_type`. A `source` column is also recommended.
#' @param blast_db Path prefix to a pre-built BLAST database used for
#'   contaminant screening.
#' @param numt_fasta Path to a FASTA file of known NUMT sequences.
#' @param run_classifier Logical. If `TRUE`, trains a Random Forest classifier
#'   on known marker types and predicts marker types for sequences labeled
#'   `"other"` or `"other_unknown"`.
#' @param run_divergence Logical. If `TRUE`, runs phylogenetic divergence
#'   checks. Requires `mafft` and `FastTree` to be installed.
#' @param run_barcode_gap Logical. If `TRUE`, runs barcode gap analysis after
#'   QC and taxonomy joining.
#' @param taxonomy_table Optional data frame used to join taxonomic lineage
#'   onto the curated sequences. Must contain `species` and the seven
#'   taxonomic ranks (kingdom through species).
#' @param ambiguity_data Optional data frame of ambiguity records to archive
#'   in the SQLite database.
#' @param db_path Optional path to write the final SQLite database.
#' @param confidence_threshold Minimum prediction probability required for
#'   the ML classifier to reclassify an unknown sequence.
#' @param min_length Minimum sequence length required for ML classification.
#' @param classifier_trees Number of trees used by the Random Forest classifier.
#' @param barcode_gap_markers Optional character vector of markers to include
#'   in barcode gap analysis. If `NULL`, all markers except `"other"` are used.
#' @param barcode_gap_min_seqs Minimum number of sequences per species
#'   required for barcode gap analysis.
#' @param barcode_gap_parallel Logical. If `TRUE`, runs barcode gap analysis
#'   in parallel where supported.
#' @param reference_table Name of the final reference table written to the
#'   SQLite database. Default `"reference_final"`.
#' @param qc_table Name of the QC table. Default `"qc_reference"`.
#' @param barcode_gap_table Name of the barcode gap table. Default `"barcode_gap_metrics"`.
#' @param ambiguity_table Name of the ambiguity table. Default `"ambiguous_sequences"`.
#' @param ambiguity_path Path to a CSV file containing ambiguity records.
#'   If `ambiguity_data = NULL`, `archive_ambiguity = TRUE`, and this file
#'   exists, it will be read and archived.
#' @param archive_ambiguity Logical. If `TRUE`, attempts to archive ambiguity
#'   records in the SQLite database.
#' @param check_ambiguity Logical. If `TRUE`, checks the input data for
#'   sequences assigned to multiple species before running QC.
#' @param on_ambiguous Action to take if ambiguous input sequences are found.
#'   One of `"stop"`, `"warn"`, or `"ignore"`.
#' @param require_taxonomy Logical. If `TRUE`, requires the final curated
#'   table to contain full kingdom-to-species taxonomy.
#' @param mafft Path to the `mafft` executable.
#' @param fasttree Path to the `FastTree` executable.
#'
#' @return An invisible list containing:
#' \describe{
#'   \item{qc1}{Data frame after first-pass QC.}
#'   \item{classifier}{Trained ML classifier object, or `NULL`.}
#'   \item{qc2}{Data frame after second-pass QC.}
#'   \item{final_data}{Final curated data frame.}
#'   \item{barcode_gap}{Barcode gap metrics, or `NULL`.}
#'   \item{db_path}{Path to the SQLite database, or `NULL`.}
#' }
#'
#' @details
#' If `check_ambiguity = TRUE`, the wrapper checks whether the same sequence
#' is assigned to multiple species before QC begins. Such ambiguities should
#' normally be resolved during standardization using `rb_resolve_ambiguous()`.
#'
#' @examples
#' \dontrun{
#' curated <- rb_curate_reference(
#'   data = standardized_data,
#'   taxonomy_table = taxonomy_table,
#'   db_path = tempfile(fileext = ".db")
#' )
#' }
#'
#' @family curation
#' @seealso `rb_run_qc_pipeline()`, `rb_build_reference_db()`
#' @export
rb_curate_reference <- function(data,
                                blast_db = NULL,
                                numt_fasta = NULL,
                                run_classifier = FALSE,
                                run_divergence = FALSE,
                                run_barcode_gap = FALSE,
                                taxonomy_table = NULL,
                                ambiguity_data = NULL,
                                db_path = NULL,
                                confidence_threshold = 0.8,
                                min_length = 100,
                                classifier_trees = 500,
                                barcode_gap_markers = NULL,
                                barcode_gap_min_seqs = 5,
                                barcode_gap_parallel = FALSE,
                                reference_table = "reference_final",
                                qc_table = "qc_reference",
                                barcode_gap_table = "barcode_gap_metrics",
                                ambiguity_table = "ambiguous_sequences",
                                ambiguity_path = "ambiguous_sequences.csv",
                                archive_ambiguity = TRUE,
                                check_ambiguity = TRUE,
                                on_ambiguous = c("stop", "warn", "ignore"),
                                require_taxonomy = TRUE,
                                mafft = "mafft",
                                fasttree = "FastTree") {
  
  rb_required_columns(data, c("species", "sequence", "seq_type"))
  
  on_ambiguous <- match.arg(on_ambiguous)
  
  if (isTRUE(check_ambiguity)){
    ambiguity_check <- rb_detect_ambiguity(data)
    n_ambiguous <- sum(ambiguity_check$is_ambiguous)
    
    if (n_ambiguous > 0) {
      msg <- sprintf("Input data contains %d ambiguous sequence(s) assigned to multiple species. Please resolve ambiguities before QC.", n_ambiguous)
      
      if (on_ambiguous == "stop") {
        stop(msg, call. = FALSE)
      }
      
      if (on_ambiguous == "warn") {
        warning(msg, call. = FALSE)
      }
    }
  }
  
  # Temporary row identifier for safe ML updates
  data$.curate_row_id <- seq_len(nrow(data))
  
  message("==================================================")
  message("Starting Reference Curation Pipeline")
  message("==================================================")
  
  # ----------------------------------------------------------
  # STAGE 1: QC Pass 1
  # ----------------------------------------------------------
  message("\n>>> STAGE 1: First-pass Quality Control")
  qc1_data <- rb_run_qc_checks(
    data,
    blast_db = blast_db,
    numt_fasta = numt_fasta,
    run_divergence = run_divergence,
    mafft = mafft,
    fasttree = fasttree
  )
  
  classifier_model <- NULL
  qc2_data <- qc1_data
  
  # ----------------------------------------------------------
  # STAGE 2: ML Classification (Optional)
  # ----------------------------------------------------------
  if (isTRUE(run_classifier)) {
    message("\n>>> STAGE 2: ML Classification of 'other' sequences")
    
    known_mask <- !qc1_data$seq_type %in% c("other", "other_unknown")
    unknown_mask <- qc1_data$seq_type %in% c("other", "other_unknown")
    
    known_table <- table(qc1_data$seq_type[known_mask])
    
    if (
      sum(known_mask) < 10 ||
      sum(unknown_mask) == 0 ||
      length(known_table) < 2 ||
      any(known_table < 2)
    ) {
      message(
        "  -> Skipping ML classification. ",
        "Need at least 10 known sequences, at least two marker classes, ",
        "and at least two sequences per class."
      )
    } else {
      message(
        sprintf(
          "  -> Training on %d known sequences across %d marker types...",
          sum(known_mask),
          length(known_table)
        )
      )
      
      known_data <- qc1_data[known_mask, ]
      
      known_features <- rb_sequence_features(known_data$sequence)
      known_features$seq_type <- factor(known_data$seq_type)
      
      tryCatch({
        classifier_model <- rb_train_classifier(
          known_features,
          label_col = "seq_type",
          trees = classifier_trees,
          mtry = 3,
          min_n = 5
        )
        
        message(sprintf("  -> Predicting %d unknown sequences...", sum(unknown_mask)))
        
        unknown_data <- qc1_data[unknown_mask, ]
        
        predictions <- rb_classify_sequences(
          classifier_model,
          unknown_data,
          confidence_threshold = confidence_threshold,
          min_length = min_length
        )
        
        high_conf_mask <- predictions$prediction_source == "ML_high_confidence"
        
        updated_rows <- predictions$.curate_row_id[high_conf_mask]
        updated_types <- predictions$final_type[high_conf_mask]
        
        found <- updated_rows %in% qc2_data$.curate_row_id
        
        if (any(found)) {
          match_idx <- match(updated_rows[found], qc2_data$.curate_row_id)
          qc2_data$seq_type[match_idx] <- updated_types[found]
        }
        
        message(sprintf("  -> Successfully reclassified %d sequences.", sum(found)))
      }, error = function(e) {
        warning("ML classification failed: ", conditionMessage(e), call. = FALSE)
      })
    }
  }
  
  # ----------------------------------------------------------
  # STAGE 3: QC Pass 2 (Re-run gene-specific checks)
  # ----------------------------------------------------------
  if (isTRUE(run_classifier) && !is.null(classifier_model)) {
    message("\n>>> STAGE 3: Second-pass Quality Control")
    qc2_data <- rb_run_qc_checks(
      qc2_data,
      blast_db = blast_db,
      numt_fasta = numt_fasta,
      run_divergence = run_divergence,
      mafft = mafft,
      fasttree = fasttree
    )
  }
  
  # ----------------------------------------------------------
  # STAGE 4: Taxonomy Joining
  # ----------------------------------------------------------
  final_data <- qc2_data
  
  if (!is.null(taxonomy_table)) {
    message("\n>>> STAGE 4: Joining Taxonomic Lineage")
    # Check taxonomy data
    if (isTRUE(require_taxonomy)) {
      rb_required_taxonomy_columns(taxonomy_table)
    }
    final_data <- rb_join_taxonomy(
      final_data,
      taxonomy_table = taxonomy_table,
      by = "species"
    )
  }
  
  if (isTRUE(require_taxonomy)) {
    rb_required_taxonomy_columns(final_data)
  }
  
  # Remove temporary row identifier before export
  qc1_data$.curate_row_id <- NULL
  qc2_data$.curate_row_id <- NULL
  final_data$.curate_row_id <- NULL
  
  barcode_gap_data <- NULL
  
  if (isTRUE(run_barcode_gap)) {
    message("\n>>> STAGE 4.5 Barcode gap analysis")
    
    if (is.null(barcode_gap_markers)) {
      barcode_gap_markers <- setdiff(unique(final_data$seq_type), c("other", "other_unknown"))
    }
    
    barcode_gap_markers <- intersect(barcode_gap_markers, unique(final_data$seq_type))
    
    if (length(barcode_gap_markers) == 0) {
      message("  -> No suitable markers available for barcode gap analysis.")
    } else {
      barcode_gap_data <- rb_run_barcode_gap(final_data, markers = barcode_gap_markers, min_seqs = barcode_gap_min_seqs, parallel = barcode_gap_parallel)
    }
  }
  
  # Read ambiguity data table
  if (is.null(ambiguity_data) && isTRUE(archive_ambiguity) && file.exists(ambiguity_path)){
    message(sprintf("  -> Reading ambiguity table from '%s'", ambiguity_path))
    ambiguity_data <- rb_read_ambiguity_csv(ambiguity_path)
  }
  
  # ----------------------------------------------------------
  # STAGE 5: SQLite Database Export
  # ----------------------------------------------------------
  if (!is.null(db_path) && nzchar(db_path)) {
    message("\n>>> STAGE 5: Writing to SQLite Database")
    rb_build_reference_db(
      path = db_path,
      reference_data = final_data,
      qc_data = qc2_data,
      barcode_gap_data = barcode_gap_data,
      ambiguity_data = ambiguity_data,
      reference_table = reference_table,
      qc_table = qc_table,
      barcode_gap_table = barcode_gap_table,
      ambiguity_table = ambiguity_table,
      require_taxonomy = require_taxonomy
    )
    message(sprintf("  -> Database saved to: %s", db_path))
  }
  
  message("\n==================================================")
  message("Curation Pipeline Complete!")
  message("==================================================")
  
  invisible(list(
    qc1 = qc1_data,
    classifier = classifier_model,
    qc2 = qc2_data,
    final_data = final_data,
    barcode_gap = barcode_gap_data,
    db_path = db_path
  ))
}

#' Run QC checks on a standardized reference table
#'
#' Internal helper used by `rb_curate_reference()`. Runs contaminant
#' screening, NUMT screening, codon validation, rRNA integrity checks,
#' genome completeness checks, optional divergence checks, and compiles
#' QC flags.
#'
#' @keywords internal
#' @noRd
rb_run_qc_checks <- function(data,
                             blast_db = NULL,
                             numt_fasta = NULL,
                             run_divergence = FALSE,
                             mafft = "mafft",
                             fasttree = "FastTree",
                             min_divergence_seqs = 5) {
  
  if (!"unique_code" %in% names(data)) {
    seq_clean <- rb_clean_sequence(data$sequence)
    unique_seqs <- unique(seq_clean)
    codes <- sprintf("seq_%05d", seq_along(unique_seqs))
    data$unique_code <- codes[match(seq_clean, unique_seqs)]
  }
  
  # 1. Contaminants
  if (!is.null(blast_db) && file.exists(paste0(blast_db, ".nsq"))) {
    message("  -> Screening contaminants...")
    data$is_contaminant <- rb_screen_contaminants(data$sequence, blast_db)
  } else {
    data$is_contaminant <- FALSE
  }
  
  # 2. NUMTs
  if (!is.null(numt_fasta) && file.exists(numt_fasta)) {
    message("  -> Screening NUMTs...")
    data$is_numt <- rb_screen_numts(data$sequence, numt_fasta)
  } else {
    data$is_numt <- FALSE
  }
  
  # 3. Codons
  message("  -> Checking codons...")
  codon_res <- mapply(
    function(seq, gene) rb_check_codons(seq, gene),
    data$sequence, data$seq_type,
    SIMPLIFY = FALSE
  )
  data$has_stop <- vapply(codon_res, function(x) ifelse(is.na(x$has_stop), FALSE, x$has_stop), logical(1))
  data$frameshifted <- vapply(codon_res, function(x) ifelse(is.na(x$frameshifted), FALSE, x$frameshifted), logical(1))
  
  # 4. rRNA integrity
  message("  -> Checking rRNA integrity...")
  rrna_res <- mapply(
    function(seq, gene) rb_check_rrna_integrity(seq, gene),
    data$sequence, data$seq_type,
    SIMPLIFY = FALSE
  )
  data$has_short <- vapply(rrna_res, function(x) ifelse(is.na(x$has_short), FALSE, x$has_short), logical(1))
  data$has_gaps <- vapply(rrna_res, function(x) ifelse(is.na(x$has_gaps), FALSE, x$has_gaps), logical(1))
  
  # 5. Genome completeness
  message("  -> Checking genome completeness...")
  data <- rb_check_genome_completeness(
    data,
    species_col = "species",
    gene_col = "seq_type"
  )
  
  # 6. Divergence
  if (isTRUE(run_divergence)) {
    message("  -> Checking divergence (requires mafft/FastTree)...")
    div_res <- tryCatch(
      rb_check_divergence(
        data,
        species_col = "species",
        gene_col = "seq_type",
        sequence_col = "sequence",
        id_col = "unique_code",
        min_seqs = min_divergence_seqs,
        mafft = mafft,
        fasttree = fasttree
      ),
      error = function(e) {
        warning("Divergence check failed: ", conditionMessage(e), call. = FALSE)
        data.frame(unique_code = data$unique_code, div_result = "processing_error", stringsAsFactors = FALSE)
      }
    )
    data <- merge(data, div_res, by = "unique_code", all.x = TRUE)
  } else {
    data$div_result <- NA_character_
  }
  
  # 7. Compile flags
  message("  -> Compiling QC flags...")
  extra_flags <- function(d) {
    vapply(seq_len(nrow(d)), function(i) {
      flags <- c()
      if (!is.na(d$genome_flag[i]) && d$genome_flag[i] == "missing_genes") flags <- c(flags, "incomplete_genome")
      if (!is.na(d$div_result[i]) && d$div_result[i] == "divergent") flags <- c(flags, "divergent")
      if (!is.na(d$div_result[i]) && d$div_result[i] %in% c("tree_error", "high_length_variation", "processing_error")) {
        flags <- c(flags, paste0("tree_", d$div_result[i]))
      }
      if (d$seq_type[i] %in% c("other", "other_unknown")) flags <- c(flags, "manual_review_needed")
      paste(flags, collapse = "|")
    }, character(1))
  }
  
  data <- rb_compile_qc_flags(data, extra_flag_fn = extra_flags)
  
  data
}
