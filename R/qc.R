# ============================================================
# R/qc.R
# Full QC pipeline: contaminant/NUMT screening, codon and rRNA
# integrity checks, genome completeness, divergence detection,
# ambiguous-content check, flag compilation, and the orchestrator that
# replaces edna_ref_qc.R AND edna_ref_qc_p2.R (call it twice — once on
# standardized data, once on ML-classified data).
# ============================================================
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
#' Screen sequences for contaminants via BLASTN against a reference db
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


# ------------------------------------------------------------
# rb_screen_numts()
                     
#' Screen sequences for NUMT contamination (query found within a
#' known NUMT reference)
#'
#' @param min_match_length Sequences shorter than this are always
#'   screened as FALSE — a short query can match a long reference by
#'   chance. Validate the false-positive rate on known-clean sequences
#'   before trusting this in production; the right threshold depends
#'   on your marker length and reference set.
# ------------------------------------------------------------

rb_screen_numts <- function(sequences, numt_fasta, min_match_length = 50) {
  if (is.null(numt_fasta) || !file.exists(numt_fasta)) {
    warning("NUMT reference file not found: ", numt_fasta, " — screening skipped.", call. = FALSE)
    return(rep(FALSE, length(sequences)))
  }
  numt_seqs <- toupper(unique(unlist(seqinr::read.fasta(numt_fasta, seqonly = TRUE))))

  # Direction: is the (short) QUERY found within a (long) NUMT
  # reference — not the reverse, since NUMT references built from
  # chromosome/scaffold/genome-scale Entrez hits are far larger than
  # any marker-length query and could never contain it the other way.
  vapply(sequences, function(seq) {
    seq_upper <- toupper(seq)
    if (nchar(seq_upper) < min_match_length) {
      return(FALSE)
    }
    any(vapply(numt_seqs, function(numt) grepl(seq_upper, numt, fixed = TRUE), logical(1)))
  }, logical(1), USE.NAMES = FALSE)
}


# ------------------------------------------------------------
# rb_check_codons()
#' Check COI/genome sequences for internal stop codons and frameshifts
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


# ------------------------------------------------------------
# rb_check_rrna_integrity()
 #' Check 12S/16S sequences for length and gap issues
# ------------------------------------------------------------

rb_check_rrna_integrity <- function(sequence, gene, rrna_genes = c("12S", "16S"),
                                     min_length = 100) {
  if (!gene %in% rrna_genes) return(list(has_short = NA, has_gaps = NA))
  list(
    has_short = nchar(sequence) < min_length,
    has_gaps = grepl("[^ACGT]", sequence)
  )
}

# ------------------------------------------------------------
# rb_check_genome_completeness()
#' Flag complete genomes missing required marker genes
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

# ------------------------------------------------------------
# rb_check_divergence()
#' Detect phylogenetically divergent (likely mislabeled) sequences
#' within each species/gene group via MAFFT + FastTree
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

#' Flag sequences with excessive ambiguous-base content (N + other
#' IUPAC codes), based on FULL sequence length
#'
#' Deliberately separate from rb_sequence_features() (R/classify.R),
#' whose GC/AT content calculation excludes ambiguous bases from its
#' denominator by design (a training-feature decision) rather than
#' penalizing for them. This is a QC-oriented check, not a feature.
rb_check_ambiguous_content <- function(sequence, max_n_fraction = 0.1) {
  seq_upper <- toupper(sequence)
  total_length <- nchar(seq_upper)
  acgt_count <- stringr::str_count(seq_upper, "[ACGT]")
  ambiguous_fraction <- ifelse(total_length > 0, 1 - (acgt_count / total_length), 0)
  ambiguous_fraction > max_n_fraction
}

# ------------------------------------------------------------
# rb_compile_qc_flags()
#' Compile per-sequence QC flags from a set of logical check columns
#' plus optional value-based extra flags
#'
#' @param checks Named list: flag name -> logical column name in `data`.
#' @param extra_flag_fn Optional `function(data) -> character vector`
#'   for flags that depend on value matching rather than a plain
#'   logical column (e.g. genome_flag == "missing_genes").
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
  flag_vec <- vapply(names(checks), function(flag_name) {
    col <- checks[[flag_name]]
    if (!col %in% names(data)) return(rep(FALSE, nrow(data)))
    val <- data[[col]]
    ifelse(is.na(val), FALSE, as.logical(val))
  }, logical(nrow(data)))

  # Unconditional reshape to matrix — handles length(checks) == 1 (where
  # vapply would otherwise return a bare vector) identically to > 1.
  flag_matrix <- matrix(flag_vec, nrow = nrow(data), ncol = length(checks),
                         dimnames = list(NULL, names(checks)))

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


# ------------------------------------------------------------
# rb_run_qc_pipeline()
           
#' Run the full QC pipeline (replaces edna_ref_qc.R AND edna_ref_qc_p2.R)
#'
#' Call this twice in user code — once on standardized data, once on
#' ML-classified data — rather than maintaining two separate scripts.
#'
#' @param parallel_workers Shared parallelism switch for contaminant
#'   screening AND the codon/rRNA integrity checks (serial when 1,
#'   furrr-parallel otherwise).
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

  if (parallel_workers > 1) {
    future::plan(future::multisession, workers = parallel_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    map2_fn <- furrr::future_map2
  } else {
    map2_fn <- purrr::map2
  }

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
  codon <- map2_fn(data[[sequence_col]], data[[gene_col]],
                    function(seq, gene) rb_check_codons(seq, gene, coding_genes = coding_genes))
  data$has_stop <- vapply(codon, function(x) isTRUE(x$has_stop), logical(1))
  data$frameshifted <- vapply(codon, function(x) isTRUE(x$frameshifted), logical(1))

  message(">>> Checking rRNA integrity")
  rrna <- map2_fn(data[[sequence_col]], data[[gene_col]],
                   function(seq, gene) rb_check_rrna_integrity(seq, gene, rrna_genes = rrna_genes))
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

