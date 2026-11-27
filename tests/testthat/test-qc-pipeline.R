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
  
  # A clean COI-like sequence (multiple of 3 (6), no internal stops in best frame)
  clean_seq <- "TCTTTATCTTTAATTGATAAGCAATGCCTTCAATCTTCAATGGATAAATATGCTGGTACTAATGGT"
  
  result <- rb_check_codons(clean_seq, gene = "COI")
  
  expect_type(result, "list")
  expect_false(is.na(result$has_stop))
  expect_false(is.na(result$frameshifted))
  expect_false(is.na(result$best_frame))
  expect_true(result$best_frame %in% c(1:6))
})


test_that("rb_check_codons does not flag partial in-frame CDS as frameshifted", {
  s <- paste0(paste(rep("ATG", 219), collapse = ""), "A")   # 658 bp, %% 3 == 1
  r <- rb_check_codons(s, gene = "COI")
  expect_false(r$frameshifted)
  expect_false(r$frame_complete)
  expect_true(r$codon_check_run)
})

test_that("genuine frameshift is detected via internal stops", {
  set.seed(42)
  orf <- synthetic_orf(400)
  s <- paste0(substr(orf, 1, 600), "A", substr(orf, 601, nchar(orf)))
  r <- rb_check_codons(s, gene = "COI")
  expect_true(r$codon_check_run)
  expect_true(r$frameshifted)
})

test_that("reverse-complemented CDS is recognised", {
  set.seed(42)
  orf <- synthetic_orf()
  rc  <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(orf)))
  expect_equal(rb_check_codons(orf, gene = "COI")$coding_strand, "F")
  expect_equal(rb_check_codons(rc,  gene = "COI")$coding_strand, "R")
})

test_that("terminal stop is not counted as internal", {
  set.seed(42)
  r <- rb_check_codons(paste0(synthetic_orf(50), "TAA"), gene = "COI")
  expect_false(r$has_stop)
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
  expect_equal(rb_compile_qc_flags(
    qc_fixture(seq_type = "16S", n_gap_chars = 1L))$qc_flag, "pass")
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

test_that("informational conditions do not block pass", {
  out <- rb_compile_qc_flags(qc_fixture(seq_type = "12S", n_ambiguous = 2L))
  expect_equal(out$qc_flag, "pass")
  expect_match(out$qc_note, "ambiguous_bases")
})

test_that("hard flags still exclude and combine", {
  out <- rb_compile_qc_flags(qc_fixture(is_contaminant = TRUE, is_numt = TRUE))
  expect_match(out$qc_flag, "contaminant")
  expect_match(out$qc_flag, "numt")
  expect_match(out$qc_flag, "\\|")
})

test_that("rb_compile_qc_flags returns 'pass' when no issues are flagged", {
  out <- rb_compile_qc_flags(qc_fixture())
  expect_equal(out$qc_flag, "pass")
  expect_true(is.na(out$qc_note))
})

test_that("rb_compile_qc_flags flags 'other' or 'other_unknown' gene as manual_review_needed", {
  out <- rb_compile_qc_flags(qc_fixture(n = 2L, seq_type = c("other", "other_unknown")))
  expect_true(all(grepl("manual_review_needed", out$qc_flag)))
})


test_that("rb_compile_qc_flags handles divergent and tree error results", {
  out <- rb_compile_qc_flags(qc_fixture(n = 3L,
                                        div_result = c("divergent", "tree_error", "pass")))
  expect_match(out$qc_flag[1], "divergent")
  expect_equal(out$qc_flag[2], "pass")            # soft → not in qc_flag
  expect_match(out$qc_note[2], "tree_error")      # and named tree_error, not tree_tree_error
  expect_equal(out$qc_flag[3], "pass")
})

test_that("sequence features use biological length", {
  r <- rb_sequence_features("ACGT-ACGTN")
  expect_equal(r$length, 9)      # 8 ACGT + 1 N; gap excluded, nchar would be 10
})

test_that("codon check runs rather than silently aborting on ambiguity", {
  set.seed(42)
  r <- rb_check_codons(paste0(synthetic_orf(30), "NNN"), gene = "COI")
  expect_true(r$codon_check_run)
})

test_that("aborted codon check is not certified as pass", {
  out <- rb_compile_qc_flags(qc_fixture(seq_type = "COI", codon_check_run = FALSE))
  expect_match(out$qc_flag, "check_failed")
})