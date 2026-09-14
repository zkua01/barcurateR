#' Assign eDNA sequences to reference taxa
#'
#' @description
#' Functions to assign amplicon sequence variants (ASVs) to a curated reference
#' database using either exact sequence matching or BLASTn alignment.
#'
#' * `rb_assign_edna()`: Main wrapper for taxonomic assignment.
#' * `rb_assign_blastn()`: Runs BLASTn against a reference database.
#' * `rb_exact_hits()`: Finds exact sequence matches in the reference.
#' * `rb_score_assignments()`: Resolves and scores assignment hits, handling ties.
#' * `rb_parse_blast_tabular()`: Parses standard BLAST tabular output.
#' * `rb_build_species_matrix()`: Aggregates ASV assignments into a species-by-sample matrix.
#' * `rb_read_fasta()`: Reads a FASTA file into a data frame (used internally for ASV input).
#'
#' @param asv_fasta Path to a FASTA file containing ASV sequences.
#' @param con A DBI connection to a barcurateR SQLite database.
#' @param db_path Path to the SQLite database (used if `con` is NULL).
#' @param marker The marker gene to extract from the database (e.g., `"12S"`, `"COI"`).
#' @param method Assignment method: `"blastn"` (uses BLAST+) or `"exact"` (dependency-free exact matching mode).
#' @param min_identity Minimum percent identity for a confident assignment.
#' @param min_coverage Minimum combined coverage (`ccovs`) for high-confidence assignment. Applied as a proportion (e.g., `0.9` for 90%). Combined coverage is calculated as the geometric mean of query and subject coverage.
#' @param max_target_seqs Maximum number of target sequences to retain per ASV in BLASTn. Defaults to 100 to ensure that ties across many congeneric species are not prematurely dropped by the BLAST algorithm.
#' @param out_dir Directory for temporary BLASTn files.
#' @param blastn Path to the `blastn` executable.
#' @param makeblastdb Path to the `makeblastdb` executable.
#' @param table_name Name of the reference table in the database (default `"reference_final"`).
#' @param hits A data frame of assignment hits (output from `rb_assign_blastn` or `rb_exact_hits`).
#' @param tie_bitscore Tolerance for bitscore ties.
#' @param asvs A data frame of ASVs (output from `rb_read_fasta`).
#' @param refs A data frame of reference sequences.
#' @param path Path to a BLAST tabular output file or FASTA file.
#' @param assignment A data frame of scored assignments.
#' @param asv_table A data frame of ASV counts per sample.
#' @param include_ambiguous Logical. Include ambiguous and `"check_full_tie"` assignments in the matrix?
#' @param unassigned Logical. Include unassigned ASVs (`no_match`) in the matrix?
#'
#' @return
#' * `rb_assign_edna()`, `rb_assign_blastn()`, `rb_exact_hits()`, and `rb_score_assignments()`
#'   return a data frame of scored assignments.
#' * `rb_parse_blast_tabular()` and `rb_read_fasta()` return a parsed data frame.
#' * `rb_build_species_matrix()` returns a species-by-sample abundance matrix.
#' 
#' @details
#' When `method = "blastn"`, BLAST+ must be installed and `blastn` and `makeblastdb` must be available on the system path.
#' 
#' The assignment pipeline utilizes advanced alignment mettrics to prevent common BLAST biases. It calculates subject coverage (`scovs`), combined coverage (`ccovs` = sqrt(qcov * scovs)), and bitscore per base (`bitscore_pb`). Hits are ranked using a robust hierarchy: `-bitscore`, `-bitscore_pb`, `evalue`, `-ccovs`, `-pident`.
#' 
#' The `assignment_status` column classifies the result:
#' \itemize{
#'   \item `"unique"`: A single best species match.
#'   \item `"ambiguous"`: Multiple species tied, but with slight metric variations.
#'   \item `"check_full_tie"`: Multiple species tied with perfectly identical top metrics across the board, requiring manual biogeographical or taxonomic curation.
#'   \item `"no_match"`: No hits passed the coverage/identity filters.
#'  }
#' 
#' For ties, the `tie_taxa` column contains a pipe-separated (`|`) string of the tied species names (e.g., `"Species A|Species B"`). If there is no tie, this column is `NA`.
#' 
#' @name rb_assignment
#' @family taxonomic assignment
NULL

#' @rdname rb_assignment
#' @export
rb_read_fasta <- function(path) {
  if (!file.exists(path)) stop("FASTA file does not exist: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  header_idx <- which(startsWith(lines, ">"))
  if (length(header_idx) == 0) stop("No FASTA records found in: ", path, call. = FALSE)
  end_idx <- c(header_idx[-1] - 1, length(lines))
  ids <- sub("^>", "", lines[header_idx])
  ids <- sub("\\s.*$", "", ids)
  seqs <- vapply(seq_along(header_idx), function(i) {
    paste(lines[(header_idx[[i]] + 1):end_idx[[i]]], collapse = "")
  }, character(1))
  data.frame(asv_id = ids, sequence = rb_clean_sequence(seqs), stringsAsFactors = FALSE)
}

#' @rdname rb_assignment
#' @export
rb_score_assignments <- function(hits, min_identity = 99, min_coverage = 0.9,
                                 tie_bitscore = 0) {
  if (nrow(hits) == 0) {
    return(data.frame())
  }
  rb_required_columns(hits, c("asv_id", "species", "pident", "qcov", "scovs", "ccovs", "evalue", "bitscore", "bitscore_pb"))
  out <- lapply(split(hits, hits$asv_id), function(x) {
    x <- x[order(-x$bitscore, -x$bitscore_pb, x$evalue, -x$ccovs, -x$pident), , drop = FALSE]
    best <- x[1, , drop = FALSE]
    tied <- x[
      x$bitscore >= best$bitscore[[1]] - tie_bitscore &
        x$bitscore_pb == best$bitscore_pb[[1]] &
        x$evalue == best$evalue[[1]] &
        x$ccovs == best$ccovs[[1]] &
        x$pident == best$pident[[1]],
      ,
      drop = FALSE
    ]
    species <- sort(unique(tied$species))
    all_metrics_identical <- length(species) > 1 &&
      nrow(tied) > 1 &&
      length(unique(paste(tied$bitscore, tied$bitscore_pb, tied$evalue, tied$ccovs, tied$pident))) == 1
    if(length(species)==1){
      status <- "unique"
    } else if(all_metrics_identical){
      status <- "check_full_tie"
    } else {
      status <- "ambiguous"
    }
    confidence <- if (status %in% c("ambiguous", "check_full_tie")) {
      "ambiguous"
    } else if (best$pident[[1]] >= min_identity && best$ccovs[[1]] >= min_coverage) {
      "high"
    } else {
      "low"
    }
    ties_out <- if(length(species)>1) paste(species, collapse = "|") else NA_character_
    data.frame(
      asv_id = best$asv_id[[1]],
      species = paste(species, collapse = "|"),
      genus = paste(sort(unique(tied$genus)), collapse = "|"),
      family = paste(sort(unique(tied$family)), collapse = "|"),
      pident = best$pident[[1]],
      qcov = best$qcov[[1]],
      scovs = best$scovs[[1]],
      ccovs = best$ccovs[[1]],
      evalue = best$evalue[[1]],
      bitscore = best$bitscore[[1]],
      bitscore_pb = best$bitscore_pb[[1]],
      n_best_species = length(species),
      assignment_status = status,
      confidence = confidence,
      tie_taxa = ties_out,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' @rdname rb_assignment
#' @export
rb_exact_hits <- function(asvs, refs) {
  ref_seq <- rb_clean_sequence(refs$sequence)
  out <- lapply(seq_len(nrow(asvs)), function(i) {
    idx <- which(ref_seq == asvs$sequence[[i]])
    if (length(idx) == 0) return(NULL)
    x <- refs[idx, , drop = FALSE]
    seq_len_val <- nchar(asvs$sequence[[i]])
    data.frame(
      asv_id = asvs$asv_id[[i]],
      species = x$species,
      genus = x$genus,
      family = x$family,
      pident = 100,
      qcov = 1,
      scovs = 1,
      ccovs = 1,
      evalue = 0,
      bitscore = seq_len_val,
      bitscore_pb = 1,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

#' @rdname rb_assignment
#' @export
rb_parse_blast_tabular <- function(path) {
  cols <- c(
    "asv_id", "subject_id", "pident", "length", "qlen", "slen",
    "evalue", "bitscore", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send"
  )
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame())
  }
  x <- utils::read.table(path, sep = "\t", header = FALSE, quote = "",
                         col.names = cols, stringsAsFactors = FALSE)
  x$qcov <- x$length / x$qlen
  x$scovs <- x$length / x$slen
  x$ccovs <- sqrt(x$qcov * x$scovs)
  x$bitscore_pb <- x$bitscore / x$length
  
  x
}

#' @rdname rb_assignment
#' @export
rb_assign_blastn <- function(asv_fasta, refs, out_dir, min_identity, min_coverage, max_target_seqs = 100, blastn, makeblastdb) {
  if (!nzchar(Sys.which(blastn))) {
    stop("BLASTN executable not found: ", blastn, ". Install BLAST+ or use method = 'exact'.", call. = FALSE)
  }
  if (!nzchar(Sys.which(makeblastdb))) {
    stop("makeblastdb executable not found: ", makeblastdb, ". Install BLAST+ or use method = 'exact'.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ref_fasta <- file.path(out_dir, "barcurateR_reference.fasta")
  db_prefix <- file.path(out_dir, "barcurateR_reference")
  blast_out <- file.path(out_dir, "barcurateR_blast.tsv")
  rb_write_fasta(paste0(">", refs$unique_code), refs$sequence, ref_fasta)
  system2(makeblastdb, c("-in", shQuote(ref_fasta), "-dbtype", "nucl", "-out", shQuote(db_prefix)),
          stdout = FALSE, stderr = FALSE)
  outfmt <- "6 qseqid sseqid pident length qlen slen evalue bitscore mismatch gapopen qstart qend sstart send"
  system2(blastn, c(
    "-query", shQuote(asv_fasta),
    "-db", shQuote(db_prefix),
    "-outfmt", shQuote(outfmt),
    "-max_target_seqs", as.character(max_target_seqs),
    "-perc_identity", as.character(min_identity),
    "-out", shQuote(blast_out)
  ), stdout = FALSE, stderr = FALSE)
  hits <- rb_parse_blast_tabular(blast_out)
  if (nrow(hits) == 0) return(data.frame())
  meta <- refs[, c("unique_code", "species", "genus", "family"), drop = FALSE]
  names(meta)[1] <- "subject_id"
  hits <- merge(hits, meta, by = "subject_id", all.x = TRUE)
  hits[hits$qcov >= min_coverage, , drop = FALSE]
}

#' @rdname rb_assignment
#' @export
rb_assign_edna <- function(asv_fasta, con = NULL, db_path = NULL, marker = NULL,
                           method = c("blastn", "exact"), min_identity = 99,
                           min_coverage = 0.9, max_target_seqs = 100,
                           out_dir = tempfile("barcurateR_assign_"),
                           blastn = "blastn", makeblastdb = "makeblastdb",
                           table_name = "reference_final") {
  method <- match.arg(method)
  own_connection <- FALSE
  if (is.null(con)) {
    con <- rb_connect(db_path)
    own_connection <- TRUE
  }
  on.exit(if (own_connection) rb_disconnect(con), add = TRUE)
  asvs <- rb_read_fasta(asv_fasta)
  refs <- rb_get_sequences(con, marker = marker, table_name = table_name)
  refs$sequence <- rb_clean_sequence(refs$sequence)
  if (method == "exact") {
    hits <- rb_exact_hits(asvs, refs)
  } else {
    hits <- rb_assign_blastn(
      asv_fasta = asv_fasta,
      refs = refs,
      out_dir = out_dir,
      min_identity = min_identity,
      min_coverage = min_coverage,
      max_target_seqs = max_target_seqs,
      blastn = blastn,
      makeblastdb = makeblastdb
    )
  }
  scored <- rb_score_assignments(hits, min_identity = min_identity, min_coverage = min_coverage)
  missing <- setdiff(asvs$asv_id, scored$asv_id)
  if (length(missing) > 0) {
    no_hits <- data.frame(
      asv_id = missing,
      species = NA_character_,
      genus = NA_character_,
      family = NA_character_,
      pident = NA_real_,
      qcov = NA_real_,
      scovs = NA_real_,
      ccovs = NA_real_,
      evalue = NA_real_,
      bitscore = NA_real_,
      bitscore_pb = NA_real_,
      n_best_species = 0L,
      assignment_status = "no_match",
      confidence = "none",
      tie_taxa = NA_character_,
      stringsAsFactors = FALSE
    )
    scored <- rbind(scored, no_hits)
  }
  scored[match(asvs$asv_id, scored$asv_id), , drop = FALSE]
}

#' @rdname rb_assignment
#' @export
rb_build_species_matrix <- function(assignment, asv_table,
                                    include_ambiguous = FALSE,
                                    unassigned = FALSE) {
  rb_required_columns(assignment, c("asv_id", "species", "assignment_status"))
  rb_required_columns(asv_table, "asv_id")
  x <- merge(assignment[, c("asv_id", "species", "assignment_status"), drop = FALSE],
             asv_table, by = "asv_id", all.x = FALSE, all.y = FALSE)
  keep <- !is.na(x$species)
  if (!include_ambiguous) keep <- keep & x$assignment_status == "unique"
  if (!unassigned) keep <- keep & x$assignment_status != "no_match"
  x <- x[keep, , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())
  sample_cols <- setdiff(names(x), c("asv_id", "species", "assignment_status"))
  for (col in sample_cols) x[[col]] <- as.numeric(x[[col]])
  stats::aggregate(x[, sample_cols, drop = FALSE], by = list(species = x$species), FUN = sum)
}
