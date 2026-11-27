test_that("off-target ambiguity is informational, not excluding", {
  s <- paste0(paste(rep("ACGT", 200), collapse = ""), "S", paste(rep("ACGT", 36), collapse = ""), "R")
  r <- rb_check_rrna_integrity(s, gene = "12S")
  expect_true(r$has_gaps)              # informational boolean
  expect_equal(r$n_ambiguous, 2L)
  expect_false(r$has_short)            # n_base large
})

test_that("length test uses biological length, not padding", {
  s <- paste0(paste(rep("A", 90), collapse = ""), paste(rep("-", 50), collapse = ""))
  r <- rb_check_rrna_integrity(s, gene = "12S", min_length = 100L)
  expect_true(r$has_short)             # n_base = 90 < 100 despite nchar = 140
})
