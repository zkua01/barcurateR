#' Screen reference sequences for contaminants and NUMTs
#'
#' @description
#' Functions to identify sequences that match known contaminants or
#' nuclear mitochondrial pseudogenes (NUMTs).
#'
#' * `rb_screen_contaminants()`: Screens sequences against a BLAST database.
#' * `rb_screen_numts()`: Screens sequences against a FASTA file of known
#'   NUMTs using exact substring matching.
#'
#' @param sequences Character vector of DNA sequences to screen.
#' @param blast_db Path prefix to a BLAST nucleotide database.
#' @param perc_identity Minimum percent identity for BLAST hits.
#' @param evalue E-value threshold for BLAST hits.
#' @param min_pident Minimum percent identity to flag as contaminant.
#' @param min_length Minimum alignment length to flag as contaminant.
#' @param threads Number of BLAST threads.
#' @param blastn Path to the `blastn` executable.
#' @param numt_fasta Path to a FASTA file of known NUMT sequences.
#' @param min_match_length Sequences shorter than this are not screened.
#'
#' @return A logical vector where `TRUE` indicates the sequence was flagged.
#' @family quality control
#' @name rb_qc_screening
NULL

#' @rdname rb_qc_screening
#' @export
rb_screen_contaminants <- function(sequences, blast_db, perc_identity = 95,
                                   evalue = 1e-20, min_pident = 95,
                                   min_length = 100, threads = 1,
                                   blastn = "blastn") {
  if (!nzchar(Sys.which(blastn))) {
    stop("blastn executable not found: ", blastn,
         ". Install BLAST+ or point `blastn` at the correct path.", call. = FALSE)
  }
  
  # Create a single temporary FASTA for all queries
  temp_fasta <- tempfile(pattern = "blast_queries_", fileext = ".fasta")
  temp_out <- tempfile(pattern = "blast_out_", fileext = ".txt")
  on.exit({
    existing <- c(temp_fasta, temp_out)[file.exists(c(temp_fasta, temp_out))]
    if (length(existing) > 0) file.remove(existing)
  }, add = TRUE)
  
  # Write all sequences to one FASTA file
  headers <- paste0(">query_", seq_along(sequences))
  lines <- character(length(sequences) * 2)
  lines[c(TRUE, FALSE)] <- headers
  lines[c(FALSE, TRUE)] <- sequences
  writeLines(lines, temp_fasta)
  
  # Run BLAST using system2() with explicit shQuote() for Windows compatibility
  blast_args <- c(
    "-query", shQuote(temp_fasta),
    "-db", shQuote(blast_db),
    "-perc_identity", as.character(perc_identity),
    "-evalue", as.character(evalue),
    "-outfmt", shQuote("6 qseqid pident evalue length"),
    "-max_target_seqs", "5",
    "-num_threads", as.character(as.integer(threads)),
    "-out", shQuote(temp_out)
  )
  
  # Capture stderr to see if BLAST complains
  blast_stderr <- system2(blastn, blast_args, stdout = FALSE, stderr = TRUE)
  
  # Initialize result vector
  is_contaminant <- rep(FALSE, length(sequences))
  
  # Parse the output file
  if (file.exists(temp_out) && file.info(temp_out)$size > 0) {
    hits <- tryCatch({
      utils::read.delim(temp_out, header = FALSE,
                 col.names = c("qseqid", "pident", "evalue", "length"),
                 colClasses = c("character", "numeric", "numeric", "integer"),
                 check.names = FALSE)
    }, error = function(e) data.frame())
    
    if (nrow(hits) > 0) {
      # Filter by user thresholds
      valid_hits <- hits[hits$pident >= min_pident & hits$length >= min_length, , drop = FALSE]
      
      # Extract the integer index from "query_1", "query_2", etc.
      indices <- as.integer(sub("^query_", "", valid_hits$qseqid))
      indices <- indices[!is.na(indices)]
      
      # Flag as contaminant
      is_contaminant[unique(indices)] <- TRUE
    }
  } else {
    # If file is empty, BLAST found no hits.
    # Only warn if BLAST produced an unexpected error in stderr 
    # (ignore the standard "0 hits found" message).
    if (length(blast_stderr) > 0 && !any(grepl("0 hits found", blast_stderr))) {
      warning("BLAST returned no hits and produced stderr:\n", 
              paste(blast_stderr, collapse = "\n"), call. = FALSE)
    }
  }
  
  is_contaminant
}

#' @rdname rb_qc_screening
#' @export
rb_screen_numts <- function(sequences, numt_fasta, min_match_length = 50) {
  if (is.null(numt_fasta) || !file.exists(numt_fasta)) {
    warning("NUMT reference file not found: ", numt_fasta, " -- screening skipped.", call. = FALSE)
    return(rep(FALSE, length(sequences)))
  }
  numt_seqs <- toupper(unique(unlist(seqinr::read.fasta(numt_fasta, seqonly = TRUE))))
  
  # Direction: is the (short) QUERY found within a (long) NUMT
  # reference -- not the reverse, since NUMT references built from
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

#' Sequence integrity checks for quality control
#'
#' @description
#' Functions to check individual sequences or groups of sequences for
#' common quality issues.
#'
#' * `rb_check_codons()`: Checks coding sequences for stop codons/frameshifts.
#' * `rb_check_rrna_integrity()`: Checks rRNA sequences for length and gaps.
#' * `rb_check_genome_completeness()`: Flags complete genomes missing markers.
#' * `rb_check_divergence()`: Detects phylogenetically divergent sequences.
#'
#' @param sequence A single DNA sequence string.
#' @param gene The gene/marker type (e.g., `"COI"`, `"12S"`).
#' @param coding_genes Gene types to check for codon issues.
#' @param genetic_code Genetic code table (default `"2"` for vertebrate mtDNA).
#' @param rrna_genes Gene types to check for rRNA integrity.
#' @param min_length Minimum expected length for rRNA sequences.
#' @param data A data frame containing sequences and metadata.
#' @param required_genes Genes required for a complete genome.
#' @param genome_gene Label for complete genome records.
#' @param species_col,gene_col,sequence_col,id_col Column names.
#' @param min_seqs Minimum sequences per group for divergence analysis.
#' @param outlier_mult Multiplier for median edge length to find outliers.
#' @param outlier_floor Minimum absolute edge length for outliers.
#' @param mafft,fasttree Paths to external executables.
#'
#' @return
#' * `rb_check_codons()` returns a list with `has_stop`, `frameshifted`, `best_frame`.
#' * `rb_check_rrna_integrity()` returns a list with `has_short`, `has_gaps`.
#' * `rb_check_genome_completeness()` returns `data` with a `genome_flag` column.
#' * `rb_check_divergence()` returns `data` with a `div_result` column.
#'
#' @family quality control
#' @name rb_qc_checks
NULL

#' @rdname rb_qc_checks
#' @export
rb_check_codons <- function(sequence, gene, coding_genes = c("COI", "complete_genome"),
                            genetic_code = "2") {
  if (!gene %in% coding_genes) {
    return(list(has_stop = NA, frameshifted = NA, best_frame = NA))
  }
  
  # Early return for empty or degenerate sequences
  seq_clean <- gsub("[^ACGT]", "N", toupper(sequence))
  if (nchar(seq_clean) == 0 || !any(grepl("[ACGT]", seq_clean))) {
    return(list(has_stop = NA, frameshifted = NA, best_frame = NA))
  }
  
  tryCatch({
    gc_table <- Biostrings::getGeneticCode(genetic_code)
    
    frame_translations <- lapply(0:2, function(frame) {
      subseq <- substr(seq_clean, frame + 1, nchar(seq_clean))
      subseq <- substr(subseq, 1, nchar(subseq) - (nchar(subseq) %% 3))
      if (nchar(subseq) == 0) return(NA)
      Biostrings::translate(Biostrings::DNAString(subseq), genetic.code = gc_table)
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

#' @rdname rb_qc_checks
#' @export
rb_check_rrna_integrity <- function(sequence, gene, rrna_genes = c("12S", "16S"),
                                    min_length = 100) {
  if (!gene %in% rrna_genes) return(list(has_short = NA, has_gaps = NA))
  list(
    has_short = nchar(sequence) < min_length,
    has_gaps = grepl("[^ACGT]", sequence)
  )
}

#' @rdname rb_qc_checks
#' @export
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

#' @rdname rb_qc_checks
#' @export
rb_check_divergence <- function(data, species_col = "species", gene_col = "gene",
                                sequence_col = "sequence", id_col = "unique_code",
                                min_seqs = 5, outlier_mult = 5, outlier_floor = 0.02,
                                mafft = "mafft", fasttree = "FastTree") {
  rb_required_columns(data, c(species_col, gene_col, sequence_col, id_col))
  
  check_one_group <- function(group_data) {
    n_seqs <- nrow(group_data)
    if (n_seqs < min_seqs) {
      return(data.frame(id = group_data[[id_col]], div_result = "insufficient_data"))
    }
    seq_lengths <- nchar(group_data[[sequence_col]])
    if (max(seq_lengths) / min(seq_lengths) > 100) {
      return(data.frame(id = group_data[[id_col]], div_result = "high_length_variation"))
    }
    missing_tools <- c(mafft, fasttree)[!nzchar(Sys.which(c(mafft, fasttree)))]
    if (length(missing_tools) > 0) {
      stop("Required external tool(s) not found: ", paste(missing_tools, collapse = ", "),
           call. = FALSE)
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

#' Compile QC flags and run the full QC pipeline
#'
#' @description
#' * `rb_compile_qc_flags()`: Combines logical check columns into a
#'   pipe-delimited `qc_flag` column.
#' * `rb_run_qc_pipeline()`: Runs the full QC pipeline in one call.
#'
#' @param data A data frame containing sequences and check results.
#' @param checks Named list mapping flag names to logical column names.
#' @param extra_flag_fn Optional function returning additional flags.
#' @param blast_db Path to a BLAST database.
#' @param numt_fasta Path to a NUMT FASTA file.
#' @param species_col,gene_col,sequence_col,id_col Column names.
#' @param min_divergence_seqs Minimum sequences for divergence analysis.
#' @param parallel_workers Number of parallel workers.
#' @param coding_genes,rrna_genes Gene types for checks.
#' @param required_genome_genes Required genes for genome completeness.
#' @param genome_gene Label for complete genome records.
#'
#' @return The input data frame with a `qc_flag` column added.
#' @family quality control
#' @name rb_qc_pipeline
NULL

#' @rdname rb_qc_pipeline
#' @export
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
  
  # Unconditional reshape to matrix -- handles length(checks) == 1 (where
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

#' @rdname rb_qc_pipeline
#' @export
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
    data[[sequence_col]], blast_db, threads = parallel_workers)

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

#' Check for excessive ambiguous-base content
#'
#' Internal QC check based on full sequence length.
#'
#' @keywords internal
#' @noRd
rb_check_ambiguous_content <- function(sequence, max_n_fraction = 0.1) {
  seq_upper <- toupper(sequence)
  total_length <- nchar(seq_upper)
  acgt_count <- stringr::str_count(seq_upper, "[ACGT]")
  ambiguous_fraction <- ifelse(total_length > 0, 1 - (acgt_count / total_length), 0)
  ambiguous_fraction > max_n_fraction
}

#' Default extra flag function for QC compilation
#'
#' Internal helper that replicates the original edna_ref_qc.R logic.
#'
#' @keywords internal
#' @noRd
default_extra_flags <- function(data) {
  vapply(seq_len(nrow(data)), function(i) {
    flags <- c()
    
    if (!is.na(data$genome_flag[i]) && data$genome_flag[i] == "missing_genes") {
      flags <- c(flags, "incomplete_genome")
    }
    if (!is.na(data$div_result[i]) && data$div_result[i] == "divergent") {
      flags <- c(flags, "divergent")
    }
    if (!is.na(data$div_result[i]) && data$div_result[i] %in% c("tree_error", "high_length_variation", "processing_error")) {
      flags <- c(flags, paste0("tree_", data$div_result[i]))
    }
    if (data$gene[i] == "other") {
      flags <- c(flags, "manual_review_needed")
    }
    
    paste(flags, collapse = "|")
  }, character(1))
}