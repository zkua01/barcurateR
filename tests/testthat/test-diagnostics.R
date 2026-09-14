# tests/testthat/test-diagnostics.R

test_that("rb_marker_coverage computes marker coverage", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_marker_coverage(con)
  expect_s3_class(res, "data.frame")
  expect_true("seq_type" %in% names(res))
  expect_true("n_sequences" %in% names(res))
  
  # Only passing sequences should be counted by default
  # COI has 2 passing, 12S has 2 passing, 16S has 1 passing
  coi_row <- res[res$seq_type == "COI", ]
  expect_equal(sum(coi_row$n_sequences), 2) 
})

test_that("rb_source_coverage computes source coverage", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_source_coverage(con)
  expect_s3_class(res, "data.frame")
  expect_true("source" %in% names(res))
})

test_that("rb_qc_summary prefers qc_reference table", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_qc_summary(con)
  expect_s3_class(res, "data.frame")
  expect_true("qc_flag" %in% names(res))
  expect_true("n_sequences" %in% names(res))
  
  # Check counts from qc_reference
  pass_count <- res$n_sequences[res$qc_flag == "pass"]
  expect_equal(sum(pass_count), 5)
})

test_that("rb_qc_summary falls back to reference_final if qc table missing", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Only write reference_final, no qc_reference
  ref_data <- data.frame(
    species = "Sp A",
    sequence = "ACGT",
    qc_flag = c("pass", "contaminant"),
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "reference_final", ref_data)
  
  res <- rb_qc_summary(con)
  expect_equal(sum(res$n_sequences), 2)
})

test_that("rb_barcode_gap reads barcode gap metrics", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_barcode_gap(con)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  
  res_filt <- rb_barcode_gap(con, species = "Species A")
  expect_equal(nrow(res_filt), 1)
})

test_that("rb_barcode_gap returns empty df if table missing", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_barcode_gap(con)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
})

test_that("rb_ambiguity reads ambiguous sequences", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_ambiguity(con)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
})