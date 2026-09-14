#' Build contaminant reference databases for QC screening
#'
#' @description
#' Functions to download contaminant sequences from NCBI and build
#' nucleotide BLAST databases for screening reference sequences against
#' known contaminants (e.g., human mtDNA, lab bacteria, NUMTs).
#'
#' * `rb_build_contaminant_db()`: Downloads sequences from NCBI based on
#'   an organisms table, saving one FASTA file per category.
#' * `rb_build_blast_db()`: Combines FASTA files and builds a nucleotide
#'   BLAST database using `makeblastdb`.
#'
#' @param organisms A data frame with columns:
#'   * `type`: Either `"accession"` (direct `entrez_fetch`) or `"query"`
#'     (`entrez_search` then fetch).
#'   * `category`: Output grouping that becomes the FASTA filename
#'     (e.g., `"human_mtDNA"`, `"lab_bacteria"`, `"numt_source"`).
#'   * `value`: The accession ID (for `type = "accession"`) or NCBI
#'     search query string (for `type = "query"`).
#' @param out_dir Directory where output FASTA files will be saved.
#'   Created automatically if it does not exist.
#' @param retmax Maximum number of sequences to retrieve per NCBI search.
#' @param requests_per_second Throttle rate between NCBI requests.
#'   The unauthenticated E-utilities limit is approximately 3 requests
#'   per second.
#' @param fasta_paths Character vector of FASTA file paths to combine
#'   and build into a BLAST database.
#' @param out_prefix Output prefix for the BLAST database files
#'   (e.g., `"data/contaminants/contaminant_db"`).
#' @param title Title string embedded in the BLAST database metadata.
#' @param makeblastdb Path to the `makeblastdb` executable.
#' @param stream_combine Logical. If `TRUE`, combines files via chunked
#'   streaming instead of loading them fully into memory. Recommended
#'   for large contaminant sets.
#'
#' @return
#' * `rb_build_contaminant_db()` returns a named character vector of
#'   output FASTA file paths, one per category.
#' * `rb_build_blast_db()` returns the `out_prefix` path invisibly.
#'
#' @details
#' The typical workflow is:
#' 1. Define an `organisms` table specifying what to download.
#' 2. Call `rb_build_contaminant_db()` to download FASTA files.
#' 3. Call `rb_build_blast_db()` on the relevant FASTA files to create
#'    the BLAST database.
#'
#' Note: NUMT sequences are typically screened via exact substring
#' matching (`rb_screen_numts()`) rather than BLAST, so NUMT FASTA files
#' are usually excluded from the BLAST database.
#'
#' @examples
#' \dontrun{
#' # Define organisms to download
#' organisms <- data.frame(
#'   type = c("accession", "query"),
#'   category = c("human_mtDNA", "lab_bacteria"),
#'   value = c(
#'     "NC_012920.1",
#'     "(Escherichia coli[ORGN] OR Pseudomonas[ORGN]) AND 16S[GENE]"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Download contaminant sequences
#' fasta_paths <- rb_build_contaminant_db(organisms, out_dir = "data/contaminants")
#'
#' # Build BLAST database (excluding NUMTs if present)
#' blast_paths <- fasta_paths[!grepl("numt", names(fasta_paths))]
#' rb_build_blast_db(blast_paths, out_prefix = "data/contaminants/contaminant_db")
#' }
#'
#' @name rb_build_contaminant_db
#' @family contaminant screening
NULL

#' @rdname rb_build_contaminant_db
#' @export
rb_build_contaminant_db <- function(organisms, out_dir, retmax = 200,
                                     requests_per_second = 3) {
  rb_required_columns(organisms, c("type", "category", "value"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  min_interval <- 1 / requests_per_second
  fetched <- list()

  for (i in seq_len(nrow(organisms))) {
    row <- organisms[i, ]
    message(sprintf("Fetching %s (%s: %s)...", row$category, row$type, row$value))
    request_start <- Sys.time()

    seqs <- if (identical(row$type, "accession")) {
      tryCatch(
        rentrez::entrez_fetch(db = "nuccore", id = row$value, rettype = "fasta"),
        error = function(e) {
          warning("Failed to fetch accession ", row$value, ": ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    } else if (identical(row$type, "query")) {
      search <- tryCatch(
        rentrez::entrez_search(db = "nuccore", term = row$value, retmax = retmax),
        error = function(e) {
          warning("Search failed for '", row$value, "': ", conditionMessage(e), call. = FALSE)
          list(ids = character(0))
        }
      )
      if (length(search$ids) == 0) {
        message("  No results found -- skipping.")
        NULL
      } else {
        tryCatch(
          rentrez::entrez_fetch(db = "nuccore", id = search$ids, rettype = "fasta"),
          error = function(e) {
            warning("Fetch failed for '", row$value, "': ", conditionMessage(e), call. = FALSE)
            NULL
          }
        )
      }
    } else {
      stop("Unknown organism `type`: '", row$type, "' (expected 'accession' or 'query')", call. = FALSE)
    }

    if (!is.null(seqs)) {
      fetched[[row$category]] <- c(fetched[[row$category]], seqs)
    }

    elapsed <- as.numeric(difftime(Sys.time(), request_start, units = "secs"))
    if (elapsed < min_interval && i < nrow(organisms)) {
      Sys.sleep(min_interval - elapsed)
    }
  }

  out_paths <- character(0)
  for (category in names(fetched)) {
    path <- file.path(out_dir, paste0(category, ".fasta"))
    writeLines(paste(fetched[[category]], collapse = "\n"), path)
    out_paths[category] <- path
    message(sprintf("  \u2713 Saved %s to '%s'", category, path))
  }
  out_paths
}
 
#' @rdname rb_build_contaminant_db
#' @export
rb_build_blast_db <- function(fasta_paths, out_prefix, title = "Contaminant Database",
                               makeblastdb = "makeblastdb", stream_combine = FALSE) {
  if (!nzchar(Sys.which(makeblastdb))) {
    stop("makeblastdb executable not found: ", makeblastdb, call. = FALSE)
  }
  missing_files <- fasta_paths[!file.exists(fasta_paths)]
  if (length(missing_files) > 0) {
    stop("Missing input FASTA file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
  }

  combined_path <- paste0(out_prefix, "_combined.fasta")
  dir.create(dirname(combined_path), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(stream_combine)) {
    out_con <- file(combined_path, open = "wb")
    on.exit(close(out_con), add = TRUE)
    for (path in fasta_paths) {
      in_con <- file(path, open = "rb")
      repeat {
        chunk <- readBin(in_con, "raw", n = 1e7)
        if (length(chunk) == 0) break
        writeBin(chunk, out_con)
      }
      close(in_con)
    }
  } else {
    contents <- unlist(lapply(fasta_paths, readLines))
    writeLines(contents, combined_path)
  }

  blast_cmd <- paste(
    makeblastdb,
    "-in", shQuote(combined_path),
    "-dbtype nucl",
    "-out", shQuote(out_prefix),
    "-title", shQuote(title)
  )
  result <- system(blast_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (result != 0) {
    stop("BLAST database creation failed. Ensure BLAST+ tools are installed and in PATH.", call. = FALSE)
  }

  db_files <- list.files(dirname(out_prefix),
                          pattern = paste0("^", basename(out_prefix), "\\..+$"),
                          full.names = TRUE)
  if (length(db_files) < 3) {
    warning("Fewer BLAST DB files created than expected. Verify database integrity.", call. = FALSE)
  }
  out_prefix
}
