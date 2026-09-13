# ============================================================
# tests/testthat/test-contaminant-screening.R
#
# Tests for contaminant database construction and screening.
#
# These tests are layered:
#
#   1. Offline validation tests.
#   2. Local BLAST+ tests, skipped if BLAST+ is not installed.
#   3. Optional live network tests, not included by default.
# ============================================================


# ------------------------------------------------------------
# rb_build_contaminant_db(): offline validation tests
# ------------------------------------------------------------

test_that("rb_build_contaminant_db requires type, category, and value columns", {
  
  bad_orgs <- data.frame(
    x = 1,
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_build_contaminant_db(bad_orgs, out_dir = tempfile()),
    "required columns"
  )
})


test_that("rb_build_contaminant_db rejects unknown organism type without needing network", {
  
  orgs <- data.frame(
    type = "bogus",
    category = "test_category",
    value = "some_value",
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_build_contaminant_db(orgs, out_dir = tempfile()),
    "Unknown organism"
  )
})


test_that("rb_build_contaminant_db handles an empty organisms table", {
  
  orgs <- data.frame(
    type = character(),
    category = character(),
    value = character(),
    stringsAsFactors = FALSE
  )
  
  out_dir <- tempfile("barcurateR_empty_contaminant_")
  
  paths <- rb_build_contaminant_db(orgs, out_dir = out_dir)
  
  expect_length(paths, 0)
})


# ------------------------------------------------------------
# rb_build_blast_db(): validation and local BLAST tests
# ------------------------------------------------------------

makeblastdb_available <- nzchar(Sys.which("makeblastdb"))


test_that("rb_build_blast_db errors when input FASTA file is missing", {
  
  skip_if_not(
    makeblastdb_available,
    "makeblastdb is not installed"
  )
  
  missing_fasta <- tempfile(fileext = ".fasta")
  
  expect_error(
    rb_build_blast_db(
      missing_fasta,
      out_prefix = tempfile("barcurateR_missing_db_")
    ),
    "Missing input FASTA"
  )
})


test_that("rb_build_blast_db errors when makeblastdb is not installed", {
  
  skip_if(
    makeblastdb_available,
    "makeblastdb is installed, so this missing-executable test is not applicable"
  )
  
  fasta_file <- tempfile(fileext = ".fasta")
  
  writeLines(
    c(
      ">test_sequence",
      "ACGTACGTACGT"
    ),
    fasta_file
  )
  
  expect_error(
    rb_build_blast_db(
      fasta_file,
      out_prefix = tempfile("barcurateR_no_makeblastdb_")
    ),
    "makeblastdb executable not found"
  )
})


test_that("rb_build_blast_db builds a tiny BLAST database from a FASTA file", {
  
  skip_if_not(
    makeblastdb_available,
    "makeblastdb is not installed"
  )
  
  skip_on_cran()
  
  tmp_dir <- tempfile("barcurateR_blast_build_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  fasta_file <- file.path(tmp_dir, "contaminant.fasta")
  
  # 120 bp sequence
  contaminant_sequence <- strrep("ACGTACGTACGT", 10)
  
  writeLines(
    c(
      ">contaminant_1",
      contaminant_sequence
    ),
    fasta_file
  )
  
  db_prefix <- file.path(tmp_dir, "contaminant_db")
  
  result <- rb_build_blast_db(
    fasta_file,
    out_prefix = db_prefix,
    title = "Test contaminant database"
  )
  
  expect_equal(result, db_prefix)
  
  db_files <- list.files(
    tmp_dir,
    pattern = "^contaminant_db\\.",
    full.names = TRUE
  )
  
  expect_gt(length(db_files), 0)
})


# ------------------------------------------------------------
# rb_screen_contaminants(): local BLAST-dependent test
# ------------------------------------------------------------

blastn_available <- nzchar(Sys.which("blastn"))


test_that("rb_screen_contaminants runs against a tiny BLAST database", {
  
  skip_if_not(
    makeblastdb_available,
    "makeblastdb is not installed"
  )
  
  skip_if_not(
    blastn_available,
    "blastn is not installed"
  )
  
  skip_on_cran()
  
  tmp_dir <- tempfile("barcurateR_screen_contaminants_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # A 120 bp contaminant sequence
  contaminant_sequence <- strrep("ACGTACGTACGT", 10)
  
  # A different 120 bp sequence
  noncontaminant_sequence <- strrep("TTTTGGGGCCCC", 10)
  
  fasta_file <- file.path(tmp_dir, "contaminant.fasta")
  
  writeLines(
    c(
      ">contaminant_1",
      contaminant_sequence
    ),
    fasta_file
  )
  
  db_prefix <- file.path(tmp_dir, "contaminant_db")
  
  rb_build_blast_db(
    fasta_file,
    out_prefix = db_prefix,
    title = "Test contaminant database"
  )
  
  query_sequences <- c(
    contaminant_sequence,
    noncontaminant_sequence
  )
  
  # At this stage, we mainly test that the function runs without error.
  # After confirming the exact return format, we can make stronger assertions.
  expect_no_error(
    result <- rb_screen_contaminants(
      query_sequences,
      blast_db = db_prefix
    )
  )
  
  expect_false(is.null(result))
})

# ------------------------------------------------------------
# rb_screen_contaminants(): local BLAST-dependent test
# ------------------------------------------------------------

test_that("rb_screen_contaminants flags exact matches and ignores non-matches", {
  
  skip_if_not(makeblastdb_available, "makeblastdb is not installed")
  skip_if_not(blastn_available, "blastn is not installed")
  skip_on_cran()
  
  tmp_dir <- tempfile("barcurateR_screen_contaminants_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Use realistic, non-repetitive 120bp sequences 
  contaminant_sequence <- "AGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAA"
  noncontaminant_sequence <- "TTAACGCCAGTTTTACCAGTCGGAGATTGGGGGTCGCGAAGTACCGTGAGGTAAAGATGCTTTGAAGAGCGTGTTATAGGGCAATAATGTCGATCTTTTATCTCTTCTT"
  
  fasta_file <- file.path(tmp_dir, "contaminant.fasta")
  writeLines(c(">contaminant_1", contaminant_sequence), fasta_file)
  
  db_prefix <- file.path(tmp_dir, "contaminant_db")
  rb_build_blast_db(fasta_file, out_prefix = db_prefix, title = "Test contaminant database")
  
  query_sequences <- c(contaminant_sequence, noncontaminant_sequence)
  
  # Run the screening
  result <- rb_screen_contaminants(
    query_sequences,
    blast_db = db_prefix,
    min_length = 50 
  )
  
  # STRICT ASSERTIONS:
  expect_type(result, "logical")
  expect_length(result, 2)
  
  # The first sequence is an exact match to the DB -> TRUE
  expect_true(result[1])
  
  # The second sequence is completely different -> FALSE
  expect_false(result[2])
})