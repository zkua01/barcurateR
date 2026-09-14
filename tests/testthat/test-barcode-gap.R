# ============================================================
# tests/testthat/test-barcode-gap.R
#
# Tests for barcode gap analysis functions:
#
#   - rb_diagnostic_sites()
#   - rb_barcode_gap_species()
#   - rb_run_barcode_gap()
#   - rb_barcode_gap()
#   - rb_read_barcode_gap()
#
# Heavier alignment tests require Biostrings and DECIPHER.
# ============================================================


# ------------------------------------------------------------
# Helper data
# ------------------------------------------------------------

make_barcode_gap_data <- function() {
  data.frame(
    species = c("Danio rerio", "Cyprinus carpio"),
    marker = c("COI", "12S"),
    gap_exists = c(TRUE, FALSE),
    gap_width = c(0.05, NA_real_),
    stringsAsFactors = FALSE
  )
}

make_barcode_gap_metrics <- function() {
  data.frame(
    species = c("Danio rerio", "Cyprinus carpio"),
    marker = c("12S", "COI"),
    n_sequences = c(10, 8),
    gap_exists = c(TRUE, FALSE),
    gap_width = c(0.045, NA_real_),
    assignment_risk = c("low", "high"),
    stringsAsFactors = FALSE
  )
}


make_aligned_species_fixture <- function() {
  data <- data.frame(
    species = c("Species A", "Species A", "Species B", "Species B"),
    genus = c("GenusA", "GenusA", "GenusB", "GenusB"),
    unique_code = c("s1", "s2", "s3", "s4"),
    stringsAsFactors = FALSE
  )
  
  aligned_seqs <- c(
    s1 = "AAAAAAAA",
    s2 = "AAAAAAAT",
    s3 = "CCCCCCCC",
    s4 = "CCCCCCCC"
  )
  
  list(
    data = data,
    aligned_seqs = aligned_seqs
  )
}


make_small_gap_input <- function() {
  a_base <- "ATGCGTACGTAGCTAGCTAGCTGATCGATC"
  b_base <- "TTGGCCAAATTTGGCCAAATTTGGCCAAAT"
  
  data.frame(
    species = rep(c("Speciesa", "Speciesb"), each = 3),
    genus = rep(c("Genusa", "Genusb"), each = 3),
    seq_type = "COI",
    sequence = c(
      a_base,
      paste0(substr(a_base, 1, 29), "A"),
      a_base,
      b_base,
      paste0(substr(b_base, 1, 29), "G"),
      b_base
    ),
    qc_flag = "pass",
    unique_code = paste0("seq_", 1:6),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# rb_diagnostic_sites()
# ------------------------------------------------------------

test_that("rb_diagnostic_sites counts fixed diagnostic positions", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  
  focal <- c("ACGT", "ACGT")
  other <- c("CCGT", "CCGT")
  
  n_diag <- rb_diagnostic_sites(focal, other)
  
  expect_equal(n_diag, 1)
})


test_that("rb_diagnostic_sites returns zero when sequences are identical", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  
  focal <- c("ACGT", "ACGT")
  other <- c("ACGT", "ACGT")
  
  n_diag <- rb_diagnostic_sites(focal, other)
  
  expect_equal(n_diag, 0)
})


# ------------------------------------------------------------
# rb_barcode_gap_species()
# ------------------------------------------------------------

test_that("rb_barcode_gap_species computes a barcode gap", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  
  fixture <- make_aligned_species_fixture()
  
  result <- rb_barcode_gap_species(
    species_name = "Species A",
    data = fixture$data,
    aligned_seqs = fixture$aligned_seqs,
    seq_len = 8,
    species_col = "species",
    genus_col = "genus",
    marker = "TEST"
  )
  
  expect_s3_class(result, "data.frame")
  
  expect_equal(result$species, "Species A")
  expect_equal(result$marker, "TEST")
  
  expect_true(result$gap_exists[[1]])
  expect_true(result$gap_width[[1]] > 0)
  
  expect_true("max_intra" %in% names(result))
  expect_true("min_inter" %in% names(result))
  expect_true("n_diagnostic_sites" %in% names(result))
  expect_true("sampling_adequacy" %in% names(result))
  expect_true("taxonomic_resolution" %in% names(result))
  expect_true("error_msg" %in% names(result))
  
  expect_true(is.na(result$error_msg[[1]]))
})


test_that("rb_barcode_gap_species returns NULL when focal species has fewer than two sequences", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  
  fixture <- make_aligned_species_fixture()
  
  expect_null(
    rb_barcode_gap_species(
      species_name = "Species C",
      data = fixture$data,
      aligned_seqs = fixture$aligned_seqs,
      seq_len = 8,
      species_col = "species",
      genus_col = "genus",
      marker = "TEST"
    )
  )
})


test_that("rb_barcode_gap_species returns NULL when there are no non-focal sequences", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  
  data_only_focal <- data.frame(
    species = c("A", "A"),
    genus = c("G", "G"),
    unique_code = c("s1", "s2"),
    stringsAsFactors = FALSE
  )
  
  aligned_only_focal <- c(
    s1 = "AAAAAAAA",
    s2 = "AAAAAAAT"
  )
  
  expect_null(
    rb_barcode_gap_species(
      species_name = "A",
      data = data_only_focal,
      aligned_seqs = aligned_only_focal,
      seq_len = 8,
      species_col = "species",
      genus_col = "genus",
      marker = "TEST"
    )
  )
})


test_that("rb_barcode_gap_species computes congeneric metrics when congeners are present", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  
  data <- data.frame(
    species = c(
      "GenusA species1",
      "GenusA species1",
      "GenusA species2",
      "GenusB species3"
    ),
    genus = c(
      "GenusA",
      "GenusA",
      "GenusA",
      "GenusB"
    ),
    unique_code = c("s1", "s2", "s3", "s4"),
    stringsAsFactors = FALSE
  )
  
  aligned_seqs <- c(
    s1 = "AAAAAAAA",
    s2 = "AAAAAAAT",
    s3 = "AAAACCCC",
    s4 = "CCCCCCCC"
  )
  
  result <- rb_barcode_gap_species(
    species_name = "GenusA species1",
    data = data,
    aligned_seqs = aligned_seqs,
    seq_len = 8,
    species_col = "species",
    genus_col = "genus",
    marker = "TEST"
  )
  
  expect_s3_class(result, "data.frame")
  
  expect_false(is.na(result$min_congeneric[[1]]))
  expect_false(is.na(result$congeneric_gap[[1]]))
  expect_true(result$n_congeners_in_db[[1]] >= 1)
})


# ------------------------------------------------------------
# rb_run_barcode_gap()
# ------------------------------------------------------------

test_that("rb_run_barcode_gap returns an empty data frame for absent markers", {
  
  input <- make_small_gap_input()
  
  result <- rb_run_barcode_gap(
    input,
    markers = "16S",
    species_col = "species",
    genus_col = "genus",
    gene_col = "seq_type",
    sequence_col = "sequence",
    qc_flag_col = "qc_flag",
    min_seqs = 2,
    parallel = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})


test_that("rb_run_barcode_gap computes metrics on a tiny dataset", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("pwalign")
  skip_on_cran()
  
  input <- make_small_gap_input()
  
  result <- rb_run_barcode_gap(
    input,
    markers = "COI",
    species_col = "species",
    genus_col = "genus",
    gene_col = "seq_type",
    sequence_col = "sequence",
    qc_flag_col = "qc_flag",
    min_seqs = 2,
    parallel = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  
  # We expect one row per species that had enough sequences.
  expect_true(nrow(result) >= 2)
  
  expect_true(all(result$marker == "COI"))
  
  expect_true(all(
    c(
      "species",
      "marker",
      "n_sequences",
      "max_intra",
      "min_inter",
      "gap_exists",
      "gap_width",
      "assignment_risk",
      "recommended_threshold",
      "error_msg"
    ) %in% names(result)
  ))
  
  # For this tiny clean dataset, we do not expect processing errors.
  expect_true(all(is.na(result$error_msg)))
  
  expect_true(
    all(result$assignment_risk %in% c("low", "moderate", "high", "failed"))
  )
})


# ------------------------------------------------------------
# rb_barcode_gap(): SQLite reading
# ------------------------------------------------------------

test_that("rb_barcode_gap returns an empty data frame when table is absent", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  result <- rb_barcode_gap(con)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})


test_that("barcode gap metrics can be written to SQLite and queried back", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  gap_metrics <- make_barcode_gap_data()
  
  rb_write_barcode_gap_table(con, gap_metrics)
  
  result <- rb_barcode_gap(con)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  
  result_coi <- rb_barcode_gap(con, marker = "COI")
  
  expect_equal(nrow(result_coi), 1)
  expect_equal(result_coi$marker, "COI")
  
  result_species <- rb_barcode_gap(con, species = "Danio rerio")
  
  expect_equal(nrow(result_species), 1)
  expect_equal(result_species$species, "Danio rerio")
})


# ------------------------------------------------------------
# rb_read_barcode_gap(): CSV reading
# ------------------------------------------------------------

test_that("rb_read_barcode_gap reads and filters a CSV file", {
  
  gap_metrics <- make_barcode_gap_metrics()
  
  tmp_csv <- tempfile(fileext = ".csv")
  utils::write.csv(gap_metrics, tmp_csv, row.names = FALSE)
  
  result <- rb_read_barcode_gap(tmp_csv)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  
  result_coi <- rb_read_barcode_gap(tmp_csv, marker = "COI")
  
  expect_equal(nrow(result_coi), 1)
  expect_equal(result_coi$marker, "COI")
  
  result_species <- rb_read_barcode_gap(tmp_csv, species = "Cyprinus carpio")
  
  expect_equal(nrow(result_species), 1)
  expect_equal(result_species$species, "Cyprinus carpio")
  
  unlink(tmp_csv)
})


test_that("rb_read_barcode_gap errors when the file does not exist", {
  
  missing_file <- tempfile(fileext = ".csv")
  
  expect_error(
    rb_read_barcode_gap(missing_file),
    "does not exist"
  )
})