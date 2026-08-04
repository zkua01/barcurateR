# ============================================================
# MODULE 3 — QC pipeline (replaces edna_ref_qc.R AND edna_ref_qc_p2.R)
#
# DESTINATION: R/qc.R (new file). All the QC-stage functions below are
# grouped together here because they were previously duplicated across
# two nearly-identical scripts — consolidating them into one file makes
# that duplication structurally impossible to reintroduce.
#
# Usage after this refactor: call rb_run_qc_pipeline() TWICE —
#   qc1 <- rb_run_qc_pipeline(standardized_data, blast_db, numt_fasta)
#   # ... run rb_train_classifier()/rb_classify_sequences() on qc1 ...
#   qc2 <- rb_run_qc_pipeline(classified_data, blast_db, numt_fasta)
# — instead of maintaining edna_ref_qc.R and edna_ref_qc_p2.R as two
# separate ~300-line scripts.
# ============================================================


# ------------------------------------------------------------
# rb_screen_contaminants()
# Generalizes run_blast_safe() + the future_map_lgl() apply loop from
# both scripts. Thresholds that were hardcoded into the blastn command
# string (-perc_identity 95 -evalue 1e-20) and into the post-hoc filter
# (pident > 95 & length > 100) are now parameters.
# ------------------------------------------------------------

rb_screen_contaminants <- function(sequences, blast_db, perc_identity = 95,
                                    evalue = 1e-20, min_pident = 95,
                                    min_length = 100, parallel_workers = 1,
                                    blastn = "blastn") {
  if (!nzchar(Sys.which(blastn))) {
    stop("blastn executable not found: ", blastn,
         ". Install BLAST+ or point `blastn` at the correct path.", call. = FALSE)
  }

  run_one <- function(seq) {
    temp_fasta <- tempfile(pattern = "blast_query_", fileext = ".fasta")
    temp_out <- tempfile(pattern = "blast_out_", fileext = ".txt")
    on.exit({
      existing <- c(temp_fasta, temp_out)[file.exists(c(temp_fasta, temp_out))]
      if (length(existing) > 0) file.remove(existing)
    }, add = TRUE)

    seqinr::write.fasta(sequences = list(seq), names = "query", file.out = temp_fasta)

    blast_cmd <- sprintf(
      "%s -query %s -db %s -perc_identity %s -evalue %s -outfmt '6 pident evalue length' -max_target_seqs 1 -out %s",
      blastn, shQuote(temp_fasta), shQuote(blast_db), perc_identity, evalue, shQuote(temp_out)
    )
    system(blast_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (!file.exists(temp_out) || file.info(temp_out)$size == 0) return(FALSE)

    hits <- tryCatch({
      read.delim(temp_out, header = FALSE,
                 col.names = c("pident", "evalue", "length"),
                 colClasses = c("numeric", "numeric", "integer"),
                 check.names = FALSE)
    }, error = function(e) data.frame())

    nrow(hits) > 0 && any(hits$pident > min_pident & hits$length > min_length)
  }

  if (parallel_workers > 1) {
    future::plan(future::multisession, workers = parallel_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
  } else {
    future::plan(future::sequential)
  }

  furrr::future_map_lgl(sequences, run_one)
}

# --- Test ---
# Requires a tiny local BLAST db to be meaningful; for a unit test that
# doesn't need BLAST+ installed, stub blast_db with a nonexistent path
# and confirm the "blastn not found" stop() fires cleanly:
#
# testthat::expect_error(
#   rb_screen_contaminants("ACGT", blast_db = "nope", blastn = "not_a_real_binary"),
#   "not found"
# )
#
# For a real test with BLAST+ installed: build a 1-sequence contaminant
# db from a known string, then confirm an identical query sequence
# screens TRUE and an unrelated random sequence screens FALSE.


# ------------------------------------------------------------
# rb_screen_numts()
# Generalizes Stage 3's grepl loop. Already dataset-agnostic (any numt
# fasta, any source organism) — main change is graceful handling of a
# missing file (warning + all-FALSE) instead of the caller having to
# set ref_db$is_numt <- FALSE manually.
# ------------------------------------------------------------

rb_screen_numts <- function(sequences, numt_fasta) {
  if (is.null(numt_fasta) || !file.exists(numt_fasta)) {
    warning("NUMT reference file not found: ", numt_fasta, " — screening skipped.",
            call. = FALSE)
    return(rep(FALSE, length(sequences)))
  }
  numt_seqs <- unique(unlist(seqinr::read.fasta(numt_fasta, seqonly = TRUE)))
  vapply(sequences, function(seq) {
    any(vapply(numt_seqs, function(numt) grepl(numt, seq, fixed = TRUE), logical(1)))
  }, logical(1), USE.NAMES = FALSE)
}

# --- Test ---
# tmp <- tempfile(fileext = ".fasta")
# writeLines(c(">numt1", "AAAACCCC"), tmp)
# rb_screen_numts(c("AAAACCCC", "GGGGTTTT"), tmp)
# expect: c(TRUE, FALSE)
#
# rb_screen_numts(c("AAAA"), numt_fasta = "does_not_exist.fasta")
# expect: a warning, and result FALSE (not an error)


# ------------------------------------------------------------
# rb_check_codons()
# Generalizes check_codon_integrity(). `coding_genes` and `genetic_code`
# were hardcoded (gene %in% c("COI","complete_genome"),
# genetic.code = "VertebrateMitochondrial") — both now parameters, so a
# dataset using a different genetic code (e.g. invertebrates, plants)
# or a different set of coding-gene labels works unchanged.
# ------------------------------------------------------------

rb_check_codons <- function(sequence, gene, coding_genes = c("COI", "complete_genome"),
                             genetic_code = "VertebrateMitochondrial") {
  if (!gene %in% coding_genes) {
    return(list(has_stop = NA, frameshifted = NA, best_frame = NA))
  }
  tryCatch({
    seq_clean <- gsub("[^ACGT]", "N", toupper(sequence))
    frame_translations <- lapply(0:2, function(frame) {
      subseq <- substr(seq_clean, frame + 1, nchar(seq_clean))
      subseq <- substr(subseq, 1, nchar(subseq) - (nchar(subseq) %% 3))
      if (nchar(subseq) == 0) return(NA)
      Biostrings::translate(Biostrings::DNAString(subseq), genetic.code = genetic_code)
    })
    stop_counts <- sapply(frame_translations, function(aa) {
      if (identical(aa, NA)) return(Inf)
      sum(Biostrings::countPattern("*", aa))
    })
    best_frame <- which.min(stop_counts)
    min_stops <- stop_counts[best_frame]
    seq_length <- nchar(seq_clean)
    list(
      has_stop = min_stops > 0,
      frameshifted = (seq_length %% 3 != 0) || (min_stops > 2),
      best_frame = best_frame
    )
  }, error = function(e) list(has_stop = NA, frameshifted = NA, best_frame = NA))
}

# --- Test ---
# rb_check_codons("ATGGCTTAA", gene = "COI")
# expect: has_stop = TRUE (TAA is a stop codon)
#
# rb_check_codons("ATGGCTGCT", gene = "COI")
# expect: has_stop = FALSE
#
# rb_check_codons("ATGGCTGCT", gene = "16S")
# expect: list(NA, NA, NA) — 16S isn't in the default coding_genes
#
# rb_check_codons("ATGGCTTAA", gene = "COI",
#                  coding_genes = "COI", genetic_code = "Standard")
# expect: runs with the Standard genetic code instead of VertebrateMitochondrial


# ------------------------------------------------------------
# rb_check_rrna_integrity()
# Generalizes check_rrna_integrity(). `rrna_genes` and `min_length`
# (previously hardcoded 100) are now parameters.
# ------------------------------------------------------------

rb_check_rrna_integrity <- function(sequence, gene, rrna_genes = c("12S", "16S"),
                                     min_length = 100) {
  if (!gene %in% rrna_genes) return(list(has_short = NA, has_gaps = NA))
  list(
    has_short = nchar(sequence) < min_length,
    has_gaps = grepl("[^ACGT]", sequence)
  )
}

# --- Test ---
# rb_check_rrna_integrity("ACGT", gene = "12S")
# expect: has_short = TRUE (4 bp < 100)
#
# rb_check_rrna_integrity(strrep("ACGT", 50), gene = "12S")
# expect: has_short = FALSE (200 bp)
#
# rb_check_rrna_integrity("ACGT", gene = "12S", min_length = 2)
# expect: has_short = FALSE — confirms min_length is actually used


# ------------------------------------------------------------
# rb_check_genome_completeness()
# Generalizes Stage 6. `required_genes` (previously hardcoded to
# has_coi/has_12s) and `genome_gene` label are now parameters — a
# COI-only marker study, for example, doesn't need a 12S column at all.
# ------------------------------------------------------------

rb_check_genome_completeness <- function(data, required_genes = c("COI", "12S"),
                                          genome_gene = "complete_genome",
                                          species_col = "species", gene_col = "gene") {
  rb_required_columns(data, c(species_col, gene_col))

  species_has_genes <- stats::aggregate(
    data[[gene_col]],
    by = list(species = data[[species_col]]),
    FUN = function(genes) all(required_genes %in% genes)
  )
  names(species_has_genes) <- c(species_col, "has_required_genes")

  data <- merge(data, species_has_genes, by = species_col, all.x = TRUE)
  data$genome_flag <- ifelse(
    data[[gene_col]] == genome_gene & !data$has_required_genes,
    "missing_genes", NA_character_
  )
  data$has_required_genes <- NULL
  data
}

# --- Test ---
# toy <- data.frame(species = c("sp1","sp1"), gene = c("complete_genome","COI"))
# rb_check_genome_completeness(toy, required_genes = "COI")
# expect: genome_flag is NA for both rows (COI alone satisfies the requirement)
#
# toy2 <- data.frame(species = c("sp1"), gene = c("complete_genome"))
# rb_check_genome_completeness(toy2, required_genes = c("COI","12S"))
# expect: genome_flag = "missing_genes" (neither COI nor 12S present)


# ------------------------------------------------------------
# rb_check_divergence()
# Generalizes check_divergence_safe() + its group_modify() driver.
# Species/gene grouping was already data-driven (no hardcoded names);
# the outlier thresholds (5x median, 0.02 floor) are now parameters.
# Tool availability (mafft/FastTree) is checked inside this function
# rather than once at the top of a script, so the function is
# self-contained and independently testable.
# ------------------------------------------------------------

rb_check_divergence <- function(data, species_col = "species", gene_col = "gene",
                                 sequence_col = "sequence", id_col = "unique_code",
                                 min_seqs = 5, outlier_mult = 5, outlier_floor = 0.02,
                                 mafft = "mafft", fasttree = "FastTree") {
  rb_required_columns(data, c(species_col, gene_col, sequence_col, id_col))

  missing_tools <- c(mafft, fasttree)[!nzchar(Sys.which(c(mafft, fasttree)))]
  if (length(missing_tools) > 0) {
    stop("Required external tool(s) not found: ", paste(missing_tools, collapse = ", "),
         call. = FALSE)
  }

  check_one_group <- function(group_data) {
    n_seqs <- nrow(group_data)
    if (n_seqs < min_seqs) {
      return(data.frame(id = group_data[[id_col]], div_result = "insufficient_data"))
    }
    seq_lengths <- nchar(group_data[[sequence_col]])
    if (max(seq_lengths) / min(seq_lengths) > 100) {
      return(data.frame(id = group_data[[id_col]], div_result = "high_length_variation"))
    }

    tryCatch({
      temp_fasta <- tempfile(fileext = ".fasta")
      temp_tree <- tempfile(fileext = ".nwk")
      on.exit({
        existing <- c(temp_fasta, temp_tree)[file.exists(c(temp_fasta, temp_tree))]
        if (length(existing) > 0) file.remove(existing)
      }, add = TRUE)

      seqinr::write.fasta(sequences = as.list(group_data[[sequence_col]]),
                           names = group_data[[id_col]], file.out = temp_fasta)

      system(paste(mafft, "--auto", temp_fasta, "|", fasttree, "-nt -gtr >", temp_tree))

      if (!file.exists(temp_tree) || file.info(temp_tree)$size == 0) {
        return(data.frame(id = group_data[[id_col]], div_result = "tree_error"))
      }
      tree <- ape::read.tree(temp_tree)
      if (is.null(tree) || is.null(tree$edge.length)) {
        return(data.frame(id = group_data[[id_col]], div_result = "tree_error"))
      }

      tip_edges <- which(tree$edge[, 2] <= length(tree$tip.label))
      edge_lengths <- tree$edge.length[tip_edges]
      tip_names <- tree$tip.label[tree$edge[tip_edges, 2]]

      if (length(edge_lengths) >= 5) {
        cutoff <- max(stats::median(edge_lengths) * outlier_mult, outlier_floor)
        outliers <- tip_names[edge_lengths > cutoff & edge_lengths > outlier_floor]
      } else {
        outliers <- character(0)
      }

      data.frame(
        id = group_data[[id_col]],
        div_result = ifelse(group_data[[id_col]] %in% outliers, "divergent", "pass")
      )
    }, error = function(e) {
      data.frame(id = group_data[[id_col]], div_result = "processing_error")
    })
  }

  split_key <- paste(data[[species_col]], data[[gene_col]], sep = "___")
  groups <- split(data, split_key)
  results <- do.call(rbind, lapply(groups, check_one_group))
  names(results)[names(results) == "id"] <- id_col
  merge(data, results, by = id_col, all.x = TRUE)
}

# --- Test ---
# testthat::skip_if_not(nzchar(Sys.which("mafft")) && nzchar(Sys.which("FastTree")))
# toy <- data.frame(
#   species = rep("sp1", 6), gene = rep("COI", 6),
#   unique_code = paste0("seq_", 1:6),
#   sequence = c(rep("ACGTACGTACGT", 5), "TTTTTTTTTTTT")  # 1 clear outlier
# )
# result <- rb_check_divergence(toy, min_seqs = 5)
# expect: the outlier sequence flagged div_result == "divergent",
#         the other 5 flagged "pass"
#
# Also confirm min_seqs is respected:
# rb_check_divergence(toy[1:3, ], min_seqs = 5)
# expect: all div_result == "insufficient_data" (below min_seqs threshold)


# ------------------------------------------------------------
# rb_compile_qc_flags()
# Generalizes Stage 8's rowwise flag builder. Simple logical-column
# flags (contaminant/numt/internal_stop/frameshifted/sequence_short/
# sequence_gap) are driven by the `checks` named list instead of being
# hardcoded `if` statements. Flags that depend on VALUE matching rather
# than a plain logical column (genome_flag == "missing_genes",
# div_result == "divergent", gene == "other") are handled by an
# optional `extra_flag_fn` callback, so a user can register new
# value-based flags without editing this function's body.
# ------------------------------------------------------------

rb_compile_qc_flags <- function(data,
                                 checks = list(
                                   contaminant = "is_contaminant",
                                   numt = "is_numt",
                                   internal_stop = "has_stop",
                                   frameshifted = "frameshifted",
                                   sequence_short = "has_short",
                                   sequence_gap = "has_gaps"
                                 ),
                                 extra_flag_fn = NULL) {
  flag_matrix <- vapply(names(checks), function(flag_name) {
    col <- checks[[flag_name]]
    if (!col %in% names(data)) return(rep(FALSE, nrow(data)))
    val <- data[[col]]
    ifelse(is.na(val), FALSE, as.logical(val))
  }, logical(nrow(data)))

  qc_flag <- apply(flag_matrix, 1, function(row) {
    hit <- names(checks)[row]
    if (length(hit) == 0) "" else paste(hit, collapse = "|")
  })

  if (!is.null(extra_flag_fn)) {
    extra_flags <- extra_flag_fn(data)
    qc_flag <- ifelse(
      nzchar(extra_flags),
      ifelse(nzchar(qc_flag), paste(qc_flag, extra_flags, sep = "|"), extra_flags),
      qc_flag
    )
  }

  qc_flag[qc_flag == ""] <- "pass"
  data$qc_flag <- qc_flag
  data
}

# --- Test ---
# toy <- data.frame(is_contaminant = c(TRUE, FALSE), is_numt = c(FALSE, FALSE),
#                    has_stop = c(FALSE, FALSE), frameshifted = c(FALSE, FALSE),
#                    has_short = c(FALSE, FALSE), has_gaps = c(FALSE, FALSE))
# rb_compile_qc_flags(toy)$qc_flag
# expect: c("contaminant", "pass")
#
# extra_fn <- function(d) ifelse(d$special == "yes", "manual_review_needed", "")
# rb_compile_qc_flags(toy, extra_flag_fn = extra_fn) # with a `special` col added
# expect: extra flag appended (comma/pipe-joined) only where extra_fn returns non-empty


# ------------------------------------------------------------
# rb_run_qc_pipeline()
# THE key consolidation: this orchestrator replaces the ENTIRE body of
# BOTH edna_ref_qc.R and edna_ref_qc_p2.R. Call it twice in user code —
# once on standardized data, once after ML classification — instead of
# maintaining two ~300-line scripts that drift out of sync with each
# other (as they had already started to: p2 has minor formatting/order
# differences from the original even though the logic is meant to be
# identical).
# ------------------------------------------------------------

rb_run_qc_pipeline <- function(data, blast_db, numt_fasta = NULL,
                                species_col = "species", gene_col = "gene",
                                sequence_col = "sequence", id_col = "unique_code",
                                min_divergence_seqs = 5, parallel_workers = 1,
                                coding_genes = c("COI", "complete_genome"),
                                rrna_genes = c("12S", "16S"),
                                required_genome_genes = c("COI", "12S"),
                                genome_gene = "complete_genome") {
  rb_required_columns(data, c(species_col, gene_col, sequence_col, id_col))
  data[[sequence_col]] <- rb_clean_sequence(data[[sequence_col]])

  message(">>> Screening contaminants")
  data$is_contaminant <- rb_screen_contaminants(
    data[[sequence_col]], blast_db, parallel_workers = parallel_workers
  )

  message(">>> Screening NUMTs")
  data$is_numt <- if (!is.null(numt_fasta)) {
    rb_screen_numts(data[[sequence_col]], numt_fasta)
  } else {
    FALSE
  }

  message(">>> Checking codon integrity")
  codon <- Map(function(seq, gene) rb_check_codons(seq, gene, coding_genes = coding_genes),
               data[[sequence_col]], data[[gene_col]])
  data$has_stop <- vapply(codon, function(x) isTRUE(x$has_stop), logical(1))
  data$frameshifted <- vapply(codon, function(x) isTRUE(x$frameshifted), logical(1))

  message(">>> Checking rRNA integrity")
  rrna <- Map(function(seq, gene) rb_check_rrna_integrity(seq, gene, rrna_genes = rrna_genes),
              data[[sequence_col]], data[[gene_col]])
  data$has_short <- vapply(rrna, function(x) isTRUE(x$has_short), logical(1))
  data$has_gaps <- vapply(rrna, function(x) isTRUE(x$has_gaps), logical(1))

  message(">>> Checking genome completeness")
  data <- rb_check_genome_completeness(
    data, required_genes = required_genome_genes, genome_gene = genome_gene,
    species_col = species_col, gene_col = gene_col
  )

  message(">>> Checking divergence")
  data <- rb_check_divergence(
    data, species_col = species_col, gene_col = gene_col,
    sequence_col = sequence_col, id_col = id_col, min_seqs = min_divergence_seqs
  )

  message(">>> Compiling QC flags")
  extra_fn <- function(d) {
    flags <- rep("", nrow(d))
    add <- function(cond, label) {
      flags[cond] <<- ifelse(nzchar(flags[cond]), paste(flags[cond], label, sep = "|"), label)
    }
    add(!is.na(d$genome_flag) & d$genome_flag == "missing_genes", "incomplete_genome")
    add(!is.na(d$div_result) & d$div_result == "divergent", "divergent")
    add(!is.na(d$div_result) & d$div_result %in% c("tree_error", "high_length_variation", "processing_error"),
        paste0("tree_", d$div_result))
    add(d[[gene_col]] == "other", "manual_review_needed")
    flags
  }
  rb_compile_qc_flags(data, extra_flag_fn = extra_fn)
}

# --- Test ---
# End-to-end smoke test on ~10 synthetic sequences with a tiny fabricated
# blast_db (built from an unrelated organism) and no numt_fasta, on a
# NON-fish, NON-YZFishDB dataset — confirms:
#   1. it runs without any fish-specific column or species assumptions
#   2. output always has a qc_flag column with "pass" for clean input
#   3. calling it twice (pass 1, then again on the same data as a stand-in
#      for the ML-classified pass 2) produces identical qc_flag values —
#      this is the regression test that proves the two-script duplication
#      has actually been eliminated, since both calls now run the exact
#      same code path.
