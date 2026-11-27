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
#' @param perc_identity Minimum percent identity passed directly to BLAST 
#'   -perc_identify flag as a hard pre-filter for BLAST hits.
#' @param evalue E-value threshold for BLAST hits.
#' @param min_pident Post-filter minimum percent identity threshold applied to 
#'   the BLAST results to flag as contaminant.
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
#' * `rb_check_codons()`: Checks coding sequences (COI) for
#'   stop codons and frameshifts.
#' * `rb_check_rrna_integrity()`: Checks rRNA sequences (12S, 16S) for abnormal
#'   length and the presence of any non-ACGT characters (including alignment gaps
#'   and ambiguous IUPAC codes like 'N').
#' * `rb_check_divergence()`: Detects phylogenetically divergent sequences
#'   within species/gene groups using MAFFT alignment and FastTree.
#'
#' @param sequence A single DNA sequence string.
#' @param gene The gene/marker type (e.g., `"COI"`, `"12S"`, `"genome"`).
#' @param coding_genes Gene types to check for codon issues.
#'   Default `c("COI")`.
#' @param genetic_code Genetic code table (default `"2"` for vertebrate mtDNA).
#' @param max_internal_stops Integer scalar (default `2L`). Maximum number of
#'   internal stop codons tolerated in the best of the six reading frames
#'   before the record is flagged `frameshifted`. The terminal stop of a
#'   complete CDS is not counted. `has_stop` is set independently whenever any
#'   internal stop is present, so lowering this value only tightens the
#'   `frameshifted` flag, not `internal_stop`.
#' @param max_ambiguous_frac Numeric scalar in `[0, 1]`, or `NULL` (default).
#'   Optional exclusion threshold on the proportion of ambiguous (IUPAC) bases
#'   among called bases (`frac_ambig`). When supplied, records with
#'   `frac_ambig > max_ambiguous_frac` receive the `ambiguous_excess` flag and
#'   are excluded from the passing set. When `NULL`, ambiguity never excludes a
#'   record; it is reported only through the informational `ambiguous_bases`
#'   note and the `n_ambiguous` / `frac_ambig` metadata columns. 
#' @param rrna_genes Gene types to check for rRNA integrity.
#' @param min_length Minimum expected length for rRNA sequences.
#' @param data A data frame containing sequences and metadata.
#' @param species_col,gene_col,sequence_col,id_col Column names.
#' @param min_seqs Minimum sequences per group for divergence analysis.
#' @param outlier_mult Multiplier for median edge length to find outliers.
#' @param outlier_floor Minimum absolute edge length for outliers.
#' @param mafft,fasttree Paths to external executables.
#'
#' @return
#' * `rb_check_codons()` returns a list with `has_stop`, `frameshifted`, `best_frame`.
#' * `rb_check_rrna_integrity()` returns a list with `has_short`, `has_gaps`.
#' * `rb_check_divergence()` returns `data` with a `div_result` column.
#' 
#' @details
#' Note: This function relies on shell pipes (|) to stream MAFFT output into FastTree. 
#' While this works on Unix-like systems, Windows users may need to ensure their shell 
#' environment (e.g., Git Bash or WSL) is properly configured, or install MSYS2, as the 
#' default cmd.exe may not support these pipes.
#' 
#'
#' @family quality control
#' @name rb_qc_checks
NULL

#' @rdname rb_qc_checks
#' @export
rb_check_codons <- function(sequence, gene,
                            coding_genes = c("COI","CYTB","complete_genome","genome","multi_marker"),
                            genetic_code = "2", max_internal_stops = 2L) {
  n <- length(sequence)
  out <- data.frame(has_stop = rep(NA, n), frameshifted = rep(NA, n),
                    n_internal_stops = rep(NA_integer_, n), best_frame = rep(NA_integer_, n),
                    coding_strand = rep(NA_character_, n), frame_complete = rep(NA, n),
                    codon_check_run = rep(NA, n), stringsAsFactors = FALSE)
  idx <- which(gene %in% coding_genes)
  if (!length(idx)) return (out)
  gc <- Biostrings::getGeneticCode(genetic_code)
  
  for(i in idx){
    out$codon_check_run[i] <- FALSE
    s <- gsub("[^ACGTRYSWKMBDHVN]", "", toupper(sequence[i]))
    if (nchar(s) < 3) next
    d <- Biostrings::DNAString(s)
    rc <- Biostrings::reverseComplement(d)
    cand <- list(F1 = d, F2 = Biostrings::subseq(d, 2), F3 = Biostrings::subseq(d, 3),
                 R1 = rc, R2 = Biostrings::subseq(rc, 2), R3 = Biostrings::subseq(rc, 3))
    cand <- cand[vapply(cand, function(x) as.integer(length(x)), integer(1)) >= 3]
    if (!length(cand)) next
    
    aa <- lapply(cand, function(x) {
      x <- Biostrings::subseq(x, 1, length(x) - (length(x) %% 3))
      tryCatch(Biostrings::translate(x, genetic.code = gc, no.init.codon = TRUE, if.fuzzy.codon = "solve"),
               error = function(e) NULL)
    })
    
    ok <- !vapply(aa, is.null, logical(1))
    if (!any(ok)) next
    aa <- aa[ok]
    nm <- names(cand)[ok]
    
    stops <- vapply(aa, function(a){
      v <- as.character(a)
      tot <- as.integer(sum(Biostrings::letterFrequency(a, "*")))
      if (nzchar(v) && substr(v, nchar(v), nchar(v)) == "*") tot - 1L else tot
    }, integer(1))
    
    b <- which.min(stops)
    out$codon_check_run[i] <- TRUE
    out$n_internal_stops[i] <- stops[b]
    out$has_stop[i] <- stops[b] > 0L
    out$best_frame[i] <- b
    out$coding_strand[i] <- substr(nm[b], 1, 1)
    out$frame_complete[i] <- (nchar(s) %% 3 == 0)
    out$frameshifted[i] <- stops[b] > max_internal_stops
  }
  out
}

#' @rdname rb_qc_checks
#' @export
rb_check_rrna_integrity <- function(sequence, gene, rrna_genes = c("12S", "16S"),
                                    min_length = 100L, max_ambiguous_frac = NULL) {
  comp <- rb_sequence_composition(sequence)
  
  n <- length(sequence)
  out <- data.frame(has_gaps = rep(NA, n),
                    has_short = rep(NA, n),
                    n_gap_chars = comp$n_dash + comp$n_dot,
                    n_ambiguous = comp$n_ambig,
                    n_base = comp$n_base,
                    frac_ambig = comp$frac_ambig,
                    max_gap_run = comp$max_gap_run, stringsAsFactors = FALSE)
  idx <- which(gene %in% rrna_genes)
  if(!length(idx)) return(out)
  out$has_gaps[idx] <- (out$n_gap_chars[idx] + comp$n_ambig[idx]) > 0
  out$has_short[idx] <- comp$n_base[idx] < min_length
  
  out
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
#' @param exclude_flags Character vector of flags that prevent `"pass"`.
#'   `NULL` (default) uses the hard-flag set defined in the function body.
#' @param max_ambiguous_frac Numeric scalar in `[0, 1]`, or `NULL` (default).
#'   Optional exclusion threshold on the proportion of ambiguous (IUPAC) bases
#'   among called bases (`frac_ambig`). When supplied, records with
#'   `frac_ambig > max_ambiguous_frac` receive the `ambiguous_excess` flag and
#'   are excluded from the passing set. When `NULL`, ambiguity never excludes a
#'   record; it is reported only through the informational `ambiguous_bases`
#'   note and the `n_ambiguous` / `frac_ambig` metadata columns. Requires the
#'   `frac_ambig` column to be present.
#' @param blast_db Path to a BLAST database.
#' @param numt_fasta Path to a NUMT FASTA file.
#' @param species_col,gene_col,sequence_col,id_col Column names.
#' @param min_divergence_seqs Minimum sequences for divergence analysis.
#' @param parallel_workers Number of parallel workers.
#' @param coding_genes,rrna_genes Gene types for checks.
#'
#' @return The input data frame with a `qc_flag` column added.
#' @family quality control
#' @name rb_qc_pipeline
NULL

#' @rdname rb_qc_pipeline
#' @export
rb_compile_qc_flags <- function(data,
                                exclude_flags = NULL,
                                max_ambiguous_frac = NULL) {
  
  hard_flags <- c("contaminant", "numt", "internal_stop", "frameshifted", "sequence_short", "divergent", "manual_review_needed", "check_failed")
  soft_flags <- c("sequence_gap", "ambiguous_bases", "frame_incomplete", "tree_error")
  
  if(is.null(exclude_flags)) exclude_flags <- hard_flags
  
  flag_map <- list(contaminant = data$is_contaminant,
                   numt = data$is_numt,
                   internal_stop = data$has_stop,
                   frameshifted = data$frameshifted,
                   sequence_short = data$has_short,
                   divergent = !is.na(data$div_result) & data$div_result == "divergent",
                   sequence_gap = data$n_gap_chars > 0,
                   ambiguous_bases = data$n_ambiguous > 0,
                   frame_incomplete = !is.na(data$frame_complete) & !data$frame_complete,
                   check_failed = !is.na(data$codon_check_run) & !data$codon_check_run)
  flag_map$tree_error <- !is.na(data$div_result) & 
    data$div_result %in% c("tree_error", "high_length_variation", "processing_error")
  flag_map$manual_review_needed <- data$seq_type %in% c("other", "other_unknown")
  if(!is.null(max_ambiguous_frac)){
    flag_map$ambiguous_excess <- !is.na(data$frac_ambig) & data$frac_ambig > max_ambiguous_frac
    exclude_flags <- union(exclude_flags, "ambiguous_excess")
  }
  
  unknown <- setdiff(names(flag_map), c(hard_flags, soft_flags, "ambiguous_excess"))
  if(length(unknown)){
    warning("Unclassified flags treated as informational: ", paste(unknown, collapse = ", "))
  }
  
  fired <- lapply(flag_map, function(x){
    x <- as.logical(x)
    x[is.na(x)] <- FALSE
    x
  })
  
  pick <- function(i, wanted){
    hits <- names(fired)[vapply(fired, function(f) f[i], logical(1))]
    hits <- intersect(hits, wanted)
    if(!length(hits)) NA_character_ else paste(hits, collapse = "|")
  }
  
  n <- nrow(data)
  data$qc_flag <- vapply(seq_len(n), function(i){
    v <- pick(i, exclude_flags)
    if(is.na(v)) "pass" else v
  }, character(1))
  data$qc_note <- vapply(seq_len(n), function(i){
    pick(i, setdiff(names(fired), exclude_flags))
  }, character(1))
  
  data
}

#' @rdname rb_qc_pipeline
#' @export
rb_run_qc_pipeline <- function(data, blast_db, numt_fasta = NULL,
                               species_col = "species", gene_col = "seq_type",
                               sequence_col = "sequence", id_col = "unique_code",
                               min_divergence_seqs = 5, parallel_workers = 1,
                               coding_genes = c("COI"),
                               rrna_genes = c("12S", "16S")) {
  rb_required_columns(data, c(species_col, gene_col, sequence_col, id_col))
  data[[sequence_col]] <- rb_clean_sequence(data[[sequence_col]])
  
  # Ensure seq_status exists
  if (!"seq_status" %in% names(data)) {
    data$seq_status <- "unknown"
  }
  
  # Remove stale QC columns from any previous pass
  stale_cols <- intersect(names(data), c(
    "div_result", "qc_flag",
    "is_contaminant", "is_numt",
    "has_stop", "frameshifted",
    "has_short", "has_gaps",
    "genome_flag"
  ))
  if (length(stale_cols) > 0) {
    data <- data[, setdiff(names(data), stale_cols), drop = FALSE]
  }
  
  if (parallel_workers > 1) {
    future::plan(future::multisession, workers = parallel_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    map2_fn <- furrr::future_map2
  } else {
    map2_fn <- purrr::map2
  }
  
  message(">>> Screening contaminants")
  data$is_contaminant <- rb_screen_contaminants(
    data[[sequence_col]], blast_db, threads = parallel_workers
  )
  
  message(">>> Screening NUMTs")
  data$is_numt <- if (!is.null(numt_fasta)) {
    rb_screen_numts(data[[sequence_col]], numt_fasta)
  } else {
    FALSE
  }
  
  message(">>> Checking codon integrity")
  codon <- map2_fn(
    data[[sequence_col]], data[[gene_col]],
    function(seq, gene) rb_check_codons(seq, gene, coding_genes = coding_genes)
  )
  data$has_stop <- vapply(codon, function(x) isTRUE(x$has_stop), logical(1))
  data$frameshifted <- vapply(codon, function(x) isTRUE(x$frameshifted), logical(1))
  
  message(">>> Checking rRNA integrity")
  rrna <- map2_fn(
    data[[sequence_col]], data[[gene_col]],
    function(seq, gene) rb_check_rrna_integrity(seq, gene, rrna_genes = rrna_genes)
  )
  data$has_short <- vapply(rrna, function(x) isTRUE(x$has_short), logical(1))
  data$has_gaps <- vapply(rrna, function(x) isTRUE(x$has_gaps), logical(1))

  message(">>> Checking divergence")
  data <- rb_check_divergence(
    data, species_col = species_col, gene_col = gene_col,
    sequence_col = sequence_col, id_col = id_col,
    min_seqs = min_divergence_seqs
  )
  
  message(">>> Compiling QC flags")
  rb_compile_qc_flags(data)
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
    if (!is.na(data$div_result[i]) && data$div_result[i] == "divergent") {
      flags <- c(flags, "divergent")
    }
    if (!is.na(data$div_result[i]) && data$div_result[i] %in% c(
      "tree_error", "high_length_variation", "processing_error"
    )) {
      flags <- c(flags, paste0("tree_", data$div_result[i]))
    }
    if ("gene" %in% names(data) && data$gene[i] == "other") {
      flags <- c(flags, "manual_review_needed")
    }
    paste(flags, collapse = "|")
  }, character(1))
}