test_that("informational conditions do not block pass", {
  d <- data.frame(
    is_contaminant = FALSE, is_numt = FALSE, has_stop = FALSE,
    frameshifted = FALSE, has_short = FALSE, div_result = "pass",
    n_gap_chars = 0L, n_ambiguous = 2L, frame_complete = NA,
    codon_check_run = NA, seq_type = "12S")
  out <- rb_compile_qc_flags(d)
  expect_equal(out$qc_flag, "pass")
  expect_match(out$qc_note, "ambiguous_bases")
})

test_that("hard flags still exclude", {
  d <- data.frame(
    is_contaminant = TRUE, is_numt = FALSE, has_stop = FALSE,
    frameshifted = FALSE, has_short = FALSE, div_result = "pass",
    n_gap_chars = 0L, n_ambiguous = 0L, frame_complete = NA,
    codon_check_run = NA, seq_type = "COI")
  expect_match(rb_compile_qc_flags(d)$qc_flag, "contaminant")
})
