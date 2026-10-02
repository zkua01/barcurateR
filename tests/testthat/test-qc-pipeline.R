# ============================================================
# tests/testthat/test-qc-pipeline.R
#
# Tests for the QC pipeline functions (stages 3–8).
# Stage 2 (contaminant screening) is tested separately in
# test-contaminant-screening.R.
# ============================================================


# ------------------------------------------------------------
# rb_screen_numts(): NUMT substring matching (Stage 3)
# ------------------------------------------------------------

test_that("rb_screen_numts detects query sequences found within NUMT references", {
  
  tmp_dir <- tempfile("numt_test_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Long NUMT reference (simulating a chromosome-scale NUMT)
  numt_seq <- paste0(
    strrep("ACGTACGTACGT", 5),   # 60 bp prefix
    "TTTTGGGGCCCCAAAA",           # unique middle section
    strrep("GGGGAAAACCCC", 5)    # 60 bp suffix
  )
  
  numt_fasta <- file.path(tmp_dir, "fish_numts.fasta")
  writeLines(c(">numt_1", numt_seq), numt_fasta)
  
  # Query that IS a substring of the NUMT (>= 50 bp)
  query_match <- strrep("ACGTACGTACGT", 5)   # 60 bp, matches prefix of numt
  
  # Query that is NOT a substring of the NUMT (>= 50 bp)
  query_no_match <- strrep("TTTTTTTTTT", 6)  # 60 bp poly-T, not in numt
  
  sequences <- c(query_match, query_no_match)
  
  result <- rb_screen_numts(sequences, numt_fasta, min_match_length = 50)
  
  expect_type(result, "logical")
  expect_length(result, 2)
  expect_true(result[1])
  expect_false(result[2])
})


test_that("rb_screen_numts returns all FALSE when NUMT file is missing", {
  
  sequences <- c("ACGTACGT", "TTTTGGGG")
  missing_file <- tempfile(fileext = ".fasta")
  
  # Should warn and return all FALSE, not error
  expect_warning(
    result <- rb_screen_numts(sequences, missing_file),
    "NUMT"
  )
  
  expect_type(result, "logical")
  expect_length(result, 2)
  expect_false(any(result))
})


test_that("rb_screen_numts returns all FALSE for empty sequences", {
  
  tmp_dir <- tempfile("numt_empty_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  numt_fasta <- file.path(tmp_dir, "fish_numts.fasta")
  writeLines(c(">numt_1", "ACGTACGT"), numt_fasta)
  
  result <- rb_screen_numts(character(0), numt_fasta)
  
  expect_type(result, "logical")
  expect_length(result, 0)
})


# ------------------------------------------------------------
# rb_check_codons(): COI codon validation (Stage 4)
# ------------------------------------------------------------

test_that("rb_check_codons returns NA for non-coding genes", {
  
  result <- rb_check_codons("ACGTACGTACGT", gene = "12S")
  
  expect_type(result, "list")
  expect_true(is.na(result$has_stop))
  expect_true(is.na(result$frameshifted))
  expect_true(is.na(result$best_frame))
})


test_that("rb_check_codons detects clean COI sequences", {
  skip_if_not_installed("Biostrings")
  
  # A clean COI-like sequence (multiple of 3, no internal stops in best frame)
  clean_seq <- "TCTTTATCTTTAATTGATAAGCAATGCCTTCAATCTTCAATGGATAAATATGCTGGTACTAATGGT"
  
  result <- rb_check_codons(clean_seq, gene = "COI")
  
  expect_type(result, "list")
  expect_false(is.na(result$has_stop))
  expect_false(is.na(result$frameshifted))
  expect_false(is.na(result$best_frame))
  expect_true(result$best_frame %in% c(1, 2, 3))
})


test_that("rb_check_codons detects frameshifted sequences", {
  
  # A sequence not divisible by 3 should be flagged as frameshifted
  frameshift_seq <- "TCTTTATCTTTAATTGATAAGCAATGCCTTCAATCTTCAATGGATAAATATGCTGGTACTAATGG"  # 65 bp (not multiple of 3)
  
  result <- rb_check_codons(frameshift_seq, gene = "COI")
  
  expect_type(result, "list")
  expect_true(result$frameshifted)
})


test_that("rb_check_codons handles empty or invalid sequences gracefully", {
  
  result <- rb_check_codons("", gene = "COI")
  expect_type(result, "list")
  expect_true(is.na(result$has_stop))
  
  result2 <- rb_check_codons("NNNNNN", gene = "COI")
  expect_type(result2, "list")
})

test_that("rb_check_codons detects stop codons in COI", {
  # COI test (unchanged)
  coi_result <- rb_check_codons("ATGAAACGCGAA", "COI")
  expect_false(coi_result$has_stop)
  
  # Non-coding markers should return NA
  rrna_result <- rb_check_codons("ATGAAACGCGAA", "12S")
  expect_true(is.na(rrna_result$has_stop))
})

# ------------------------------------------------------------
# rb_check_rrna_integrity(): 12S/16S checks (Stage 5)
# ------------------------------------------------------------

test_that("rb_check_rrna_integrity returns NA for non-rRNA genes", {
  
  result <- rb_check_rrna_integrity("ACGTACGTACGT", gene = "COI")
  
  expect_type(result, "list")
  expect_true(is.na(result$has_short))
  expect_true(is.na(result$has_gaps))
})


test_that("rb_check_rrna_integrity flags short sequences", {
  
  short_seq <- strrep("A", 50)  # 50 bp, below 100 threshold
  result <- rb_check_rrna_integrity(short_seq, gene = "12S")
  
  expect_type(result, "list")
  expect_true(result$has_short)
  expect_false(result$has_gaps)
})


test_that("rb_check_rrna_integrity flags sequences with gaps or ambiguity", {
  
  gapped_seq <- strrep("A", 150)
  gapped_seq <- paste0(substr(gapped_seq, 1, 75), "-", substr(gapped_seq, 76, 150))
  
  result <- rb_check_rrna_integrity(gapped_seq, gene = "16S")
  
  expect_type(result, "list")
  expect_false(result$has_short)
  expect_true(result$has_gaps)
})


test_that("rb_check_rrna_integrity passes clean sequences", {
  
  clean_seq <- strrep("ACGT", 50)  # 200 bp, no gaps
  result <- rb_check_rrna_integrity(clean_seq, gene = "12S")
  
  expect_type(result, "list")
  expect_false(result$has_short)
  expect_false(result$has_gaps)
})


# ------------------------------------------------------------
# rb_compile_qc_flags(): flag compilation (Stage 8)
# ------------------------------------------------------------

test_that("rb_compile_qc_flags returns 'pass' when no issues are flagged", {
  
  test_data <- data.frame(
    is_contaminant = FALSE,
    is_numt = FALSE,
    has_stop = FALSE,
    frameshifted = FALSE,
    has_short = FALSE,
    has_gaps = FALSE,
    genome_flag = NA_character_,
    div_result = NA_character_,
    gene = "COI",
    stringsAsFactors = FALSE
  )
  
  result <- rb_compile_qc_flags(test_data)
  
  expect_s3_class(result, "data.frame")
  expect_true("qc_flag" %in% names(result))
  expect_equal(result$qc_flag, "pass")
})


test_that("rb_compile_qc_flags combines multiple flags with pipe separator", {
  
  test_data <- data.frame(
    is_contaminant = TRUE,
    is_numt = TRUE,
    has_stop = FALSE,
    frameshifted = FALSE,
    has_short = FALSE,
    has_gaps = FALSE,
    genome_flag = NA_character_,
    div_result = NA_character_,
    gene = "COI",
    stringsAsFactors = FALSE
  )
  
  result <- rb_compile_qc_flags(test_data)
  
  expect_true(grepl("contaminant", result$qc_flag))
  expect_true(grepl("numt", result$qc_flag))
  expect_true(grepl("\\|", result$qc_flag))
})


test_that("rb_compile_qc_flags flags 'other' gene as manual_review_needed", {
  
  test_data <- data.frame(
    is_contaminant = FALSE,
    is_numt = FALSE,
    has_stop = FALSE,
    frameshifted = FALSE,
    has_short = FALSE,
    has_gaps = FALSE,
    genome_flag = NA_character_,
    div_result = NA_character_,
    gene = "other",
    stringsAsFactors = FALSE
  )
  
  result <- rb_compile_qc_flags(test_data, extra_flag_fn = default_extra_flags)
  
  expect_true(grepl("manual_review_needed", result$qc_flag))
})


test_that("rb_compile_qc_flags handles divergent and tree error results", {
  
  test_data <- data.frame(
    is_contaminant = FALSE,
    is_numt = FALSE,
    has_stop = FALSE,
    frameshifted = FALSE,
    has_short = FALSE,
    has_gaps = FALSE,
    genome_flag = NA_character_,
    div_result = c("divergent", "tree_error", "pass"),
    gene = c("COI", "COI", "COI"),
    stringsAsFactors = FALSE
  )
  
  result <- rb_compile_qc_flags(test_data, extra_flag_fn = default_extra_flags)
  
  expect_true(grepl("divergent", result$qc_flag[1]))
  expect_true(grepl("tree_tree_error", result$qc_flag[2]))
  expect_equal(result$qc_flag[3], "pass")
})
