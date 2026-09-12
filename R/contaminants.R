# ============================================================
# R/contaminants.R
# Contaminant reference database building (NCBI fetch + BLAST db build).
# ============================================================
 
# ------------------------------------------------------------
# rb_build_contaminant_db()

#' Build contaminant reference FASTAs from an organisms table
#'
#' @param organisms Data frame with columns:
#'   - type: "accession" (direct entrez_fetch) or "query" (entrez_search
#'     then fetch)
#'   - category: output grouping (e.g. "reference_contaminant",
#'     "lab_contaminant", "numt_source") — becomes the output filename
#'   - value: the accession ID (type = "accession") or NCBI search
#'     query string (type = "query")
#' @param requests_per_second Fixed-interval throttle between NCBI
#'   requests (unauthenticated E-utilities cap is ~3/sec).
#' @return Named character vector of output FASTA paths, one per
#'   category.

# Generalizes all 3 hardcoded blocks (human mtDNA via a fixed
# accession, fish NUMT source organisms via a fixed Entrez query, lab
# bacteria via another fixed query) into one function driven by an
# `organisms` table. Each row is either a direct accession fetch or an
# Entrez search+fetch, tagged with a `category` that determines which
# output FASTA it lands in — so the original's 3 fixed categories
# (human mtDNA / fish NUMTs / lab bacteria) become just 3 example rows
# in a caller-supplied table, not 3 hardcoded code paths.
# ------------------------------------------------------------
 
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
        message("  No results found — skipping.")
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
 
 
 
# ------------------------------------------------------------
# rb_build_blast_db()
#' Combine FASTA files and build a nucleotide BLAST database
#'
#' @param fasta_paths Character vector of FASTA file paths to combine.
#' @param stream_combine If TRUE, combines files via chunked
#'   file-to-file streaming instead of loading them fully into R
#'   memory — recommended for large/user-supplied contaminant sets.

# Generalizes Stage 4's makeblastdb system call. `fasta_paths` is now a
# list instead of the 2 hardcoded files (human_mtDNA.fasta +
# lab_bacteria.fasta) — the caller decides which of
# rb_build_contaminant_db()'s outputs to include (e.g. deliberately
# excluding numt_source.fasta, matching the original script's intent
# that NUMTs are screened separately via rb_screen_numts(), not BLAST).
# ------------------------------------------------------------
 
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
