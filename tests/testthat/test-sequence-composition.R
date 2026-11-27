test_that("composition separates padding, ambiguity and bases", {
  x <- rb_sequence_composition(c("ACGTACGT", "ACGT-ACGT", "ACGTNACGT", "ACGT..ACGT"))
  expect_equal(x$n_base, c(8L, 8L, 9L, 8L))          # gaps excluded from n_base
  expect_equal(x$n_dash, c(0L, 1L, 0L, 0L))
  expect_equal(x$n_dot,  c(0L, 0L, 0L, 2L))
  expect_equal(x$n_ambig, c(0L, 0L, 1L, 0L))
  expect_equal(x$n_char, c(8L, 9L, 9L, 10L))
  expect_true(is.na(rb_sequence_composition("")$frac_ambig))
})
