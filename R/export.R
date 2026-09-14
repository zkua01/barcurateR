#' Export curated reference sequences to common metabarcoding formats
#'
#' @description
#' Functions to export curated reference sequences and taxonomy to
#' FASTA-based formats used by common metabarcoding pipelines.
#'
#' * `rb_export_fasta()`: General dispatcher that writes to one of the
#'   supported formats.
#' * `rb_export_blastn()`: Writes BLAST-style FASTA headers.
#' * `rb_export_dada2()`: Writes DADA2-style FASTA headers.
#' * `rb_export_qiime2()`: Writes a QIIME2-compatible FASTA file and a
#'   separate taxonomy TSV file.
#' * `rb_export_usearch()`: Writes USEARCH/VSEARCH-style FASTA headers.
#'
#' @param data A data frame containing at least `sequence` and the seven
#'   taxonomic ranks: `kingdom`, `phylum`, `class`, `order`, `family`,
#'   `genus`, and `species`.
#' @param file Output file path. Used by `rb_export_fasta()`.
#' @param format Export format. One of `"blastn"`, `"dada2"`, `"qiime2"`,
#'   or `"usearch"`. Used by `rb_export_fasta()`.
#' @param fasta_file Output FASTA file path. Used by the format-specific
#'   export functions.
#' @param taxonomy_file Output taxonomy TSV file path. Required for
#'   `rb_export_qiime2()`.
#'
#' @return
#' * `rb_export_fasta()`, `rb_export_blastn()`, `rb_export_dada2()`, and
#'   `rb_export_usearch()` return the output FASTA path invisibly.
#' * `rb_export_qiime2()` returns an invisible list containing `fasta`
#'   and `taxonomy` file paths.
#'
#' @details
#' All export functions require full taxonomy from kingdom to species.
#' They do not perform QC filtering themselves.
#'
#' To export only QC-passing sequences, filter the data first using
#' [rb_get_sequences()] or [rb_filter_reference()]. For example:
#'
#' ```r
#' pass_refs <- rb_get_sequences(con, qc_flag = "pass")
#' rb_export_dada2(pass_refs, "reference_dada2.fasta")
#' ```
#'
#' Header styles differ by format:
#'
#' | Format | Example header style |
#' |---|---|
#' | BLAST | `>seq_id; Kingdom; Phylum; Class; ...` |
#' | DADA2 | `>seq_id;tax=Kingdom;Phylum;Class;...` |
#' | USEARCH | `>seq_id tax=Kingdom;Phylum;Class;...` |
#' | QIIME2 | Plain FASTA headers plus separate taxonomy TSV |
#'
#' @examples
#' \dontrun{
#' con <- rb_connect()
#' refs <- rb_get_sequences(con, qc_flag = "pass")
#' rb_disconnect(con)
#'
#' rb_export_dada2(refs, "ref_dada2.fasta")
#' rb_export_qiime2(refs, "ref_qiime2.fasta", "ref_taxonomy.tsv")
#' }
#'
#' @name rb_export
#' @family export
NULL

#' Write FASTA lines to file
#'
#' Internal helper used by the export functions.
#'
#' @keywords internal
#' @noRd
rb_write_fasta <- function(headers, sequences, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  lines <- character(length(headers) * 2)
  lines[seq(1, length(lines), by = 2)] <- headers
  lines[seq(2, length(lines), by = 2)] <- rb_clean_sequence(sequences)
  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}

#' Generate syntactically valid sequence identifiers
#'
#' Internal helper used by the export functions.
#'
#' @keywords internal
#' @noRd
rb_safe_id <- function(data) {
  if ("unique_code" %in% names(data)) {
    make.names(data$unique_code, unique = TRUE)
  } else if ("sequence_id" %in% names(data)) {
    make.names(data$sequence_id, unique = TRUE)
  } else {
    paste0("seq_", seq_len(nrow(data)))
  }
}

#' @rdname rb_export
#' @export
rb_export_fasta <- function(data, file, format = c("blastn", "dada2", "qiime2", "usearch"),
                            taxonomy_file = NULL) {
  format <- match.arg(format)
  switch(
    format,
    blastn = rb_export_blastn(data, file),
    dada2 = rb_export_dada2(data, file),
    qiime2 = rb_export_qiime2(data, file, taxonomy_file),
    usearch = rb_export_usearch(data, file)
  )
}

#' @rdname rb_export
#' @export
rb_export_blastn <- function(data, fasta_file) {
  rb_required_columns(data, c("sequence", "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  tax <- rb_build_taxonomy_string(data, style = "plain")
  headers <- paste0(">", rb_safe_id(data), ";", tax)
  rb_write_fasta(headers, data$sequence, fasta_file)
}

#' @rdname rb_export
#' @export
rb_export_dada2 <- function(data, fasta_file) {
  rb_required_columns(data, c("sequence", "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  tax <- rb_build_taxonomy_string(data, style = "plain")
  headers <- paste0(">", rb_safe_id(data), ";tax=", tax)
  rb_write_fasta(headers, data$sequence, fasta_file)
}

#' @rdname rb_export
#' @export
rb_export_qiime2 <- function(data, fasta_file, taxonomy_file) {
  if (is.null(taxonomy_file)) stop("taxonomy_file is required for QIIME2 export.", call. = FALSE)
  rb_required_columns(data, c("sequence", "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  ids <- rb_safe_id(data)
  rb_write_fasta(paste0(">", ids), data$sequence, fasta_file)
  tax <- rb_build_taxonomy_string(data, style = "rank_prefix")
  lines <- c("Feature ID\tTaxon", paste(ids, tax, sep = "\t"))
  dir.create(dirname(taxonomy_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, taxonomy_file, useBytes = TRUE)
  invisible(list(fasta = fasta_file, taxonomy = taxonomy_file))
}

#' @rdname rb_export
#' @export
rb_export_usearch <- function(data, fasta_file) {
  rb_required_columns(data, c("sequence", "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  tax <- rb_build_taxonomy_string(data, style = "rank_prefix")
  headers <- paste0(">", rb_safe_id(data), " tax=", tax)
  rb_write_fasta(headers, data$sequence, fasta_file)
}
