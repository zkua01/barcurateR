
 
# ------------------------------------------------------------
# rb_build_contaminant_db()
# Generalizes all 3 hardcoded blocks (human mtDNA via a fixed
# accession, fish NUMT source organisms via a fixed Entrez query, lab
# bacteria via another fixed query) into one function driven by an
# `organisms` table. Each row is either a direct accession fetch or an
# Entrez search+fetch, tagged with a `category` that determines which
# output FASTA it lands in — so the original's 3 fixed categories
# (human mtDNA / fish NUMTs / lab bacteria) become just 3 example rows
# in a caller-supplied table, not 3 hardcoded code paths.
# ------------------------------------------------------------
 
rb_build_contaminant_db <- function(organisms, out_dir, retmax = 200) {
  rb_required_columns(organisms, c("type", "category", "value"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
 
  fetched <- list()
 
  for (i in seq_len(nrow(organisms))) {
    row <- organisms[i, ]
    message(sprintf("Fetching %s (%s: %s)...", row$category, row$type, row$value))
 
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
 
# --- Test ---
# NOTE: requires network access + rentrez, and can't be fully offline-tested
# without mocking entrez_fetch/entrez_search. For CI, wrap these calls
# behind an injectable `fetcher`/`searcher` argument if you want a fully
# mockable unit test — flag if you'd like that added.
#
# Structural test (offline, using an organisms table with deliberately
# NON-fish, NON-YZFishDB categories):
#
# orgs <- data.frame(
#   type     = c("accession", "query"),
#   category = c("reference_contaminant", "lab_contaminant"),
#   value    = c("NC_012920.1", "Escherichia coli[ORGN] AND 16S[GENE]"),
#   stringsAsFactors = FALSE
# )
# paths <- rb_build_contaminant_db(orgs, out_dir = tempdir())
# expect: paths named "reference_contaminant" and "lab_contaminant"
#         (not "human_mtDNA"/"lab_bacteria" — proves category naming
#         isn't hardcoded), each pointing at a real fasta file
#
# Also test the unknown-type error path (no network needed):
# testthat::expect_error(
#   rb_build_contaminant_db(data.frame(type = "bogus", category = "x", value = "y"),
#                            out_dir = tempdir()),
#   "Unknown organism"
# )
 
 
# ------------------------------------------------------------
# rb_build_blast_db()
# Generalizes Stage 4's makeblastdb system call. `fasta_paths` is now a
# list instead of the 2 hardcoded files (human_mtDNA.fasta +
# lab_bacteria.fasta) — the caller decides which of
# rb_build_contaminant_db()'s outputs to include (e.g. deliberately
# excluding numt_source.fasta, matching the original script's intent
# that NUMTs are screened separately via rb_screen_numts(), not BLAST).
# ------------------------------------------------------------
 
rb_build_blast_db <- function(fasta_paths, out_prefix, title = "Contaminant Database",
                               makeblastdb = "makeblastdb") {
  if (!nzchar(Sys.which(makeblastdb))) {
    stop("makeblastdb executable not found: ", makeblastdb, call. = FALSE)
  }
  missing_files <- fasta_paths[!file.exists(fasta_paths)]
  if (length(missing_files) > 0) {
    stop("Missing input FASTA file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
  }
 
  combined_path <- paste0(out_prefix, "_combined.fasta")
  dir.create(dirname(combined_path), recursive = TRUE, showWarnings = FALSE)
 
  contents <- unlist(lapply(fasta_paths, readLines))
  writeLines(contents, combined_path)
 
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
 
# --- Test ---
# testthat::skip_if_not(nzchar(Sys.which("makeblastdb")))
# f1 <- tempfile(fileext = ".fasta"); writeLines(c(">a", "ACGTACGT"), f1)
# f2 <- tempfile(fileext = ".fasta"); writeLines(c(">b", "TTTTGGGG"), f2)
# prefix <- file.path(tempdir(), "test_db")
# rb_build_blast_db(list(f1, f2), out_prefix = prefix, title = "Test DB")
# expect: prefix.nhr / .nin / .nsq (or newer BLAST+ equivalents) created
#
# testthat::expect_error(
#   rb_build_blast_db("does_not_exist.fasta", out_prefix = prefix),
#   "Missing input FASTA"
# )
