# ============================================================
# tests/testthat/test-curate.R
#
# Tests for the high-level rb_curate_reference() wrapper.
#
# Updated to include:
#   - taxonomy table (required)
#   - ambiguity checking before QC
#   - optional barcode gap analysis
#   - ambiguity archiving from CSV
#   - generic table names
# ============================================================


# ------------------------------------------------------------
# Helper: test data
# ------------------------------------------------------------

make_wrapper_test_data <- function() {
  
  set.seed(42)
  
  # Track all generated sequences so that no two species share the same sequence.
  seen <- character(0)
  
  make_unique_sequences <- function(n, length) {
    out <- character(n)
    
    for (i in seq_len(n)) {
      repeat {
        candidate <- paste(
          sample(c("A", "C", "G", "T"), length, replace = TRUE),
          collapse = ""
        )
        
        if (!candidate %in% seen) break
      }
      
      out[i] <- candidate
      seen <<- c(seen, candidate)
    }
    
    out
  }
  
  # COI-like sequences: length 300
  coi_seqs <- make_unique_sequences(10, 300)
  
  # 12S-like sequences: length 150
  s12_seqs <- make_unique_sequences(10, 150)
  
  # Unknown sequences that resemble COI in length but are still unique
  other_seqs <- make_unique_sequences(5, 300)
  
  data.frame(
    species = paste0("Species_", 1:25),
    sequence = c(coi_seqs, s12_seqs, other_seqs),
    seq_type = c(
      rep("COI", 10),
      rep("12S", 10),
      rep("other", 5)
    ),
    source = "test",
    stringsAsFactors = FALSE
  )
}

make_wrapper_taxonomy <- function(species) {
  data.frame(
    species = species,
    genus = paste0("Genus_", seq_along(species)),
    family = "Testidae",
    order = "Testorder",
    class = "Testclass",
    phylum = "Testphylum",
    kingdom = "Testkingdom",
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Basic wrapper: no ML, no barcode gap
# ------------------------------------------------------------

test_that("rb_curate_reference runs QC and exports SQLite without ML", {
  
  skip_if_not_installed("RSQLite")
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  db_path <- tempfile(fileext = ".db")
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE,
    db_path = db_path
  )
  
  expect_type(result, "list")
  expect_true(all(c("qc1", "qc2", "final_data", "db_path") %in% names(result)))
  
  # QC flags should be compiled
  expect_true("qc_flag" %in% names(result$qc1))
  
  # Temporary row ID should be removed
  expect_false(".curate_row_id" %in% names(result$final_data))
  expect_false(".curate_row_id" %in% names(result$qc1))
  expect_false(".curate_row_id" %in% names(result$qc2))
  
  # Taxonomy columns should be present
  expect_true("kingdom" %in% names(result$final_data))
  expect_true("genus" %in% names(result$final_data))
  
  # SQLite database should exist and contain the generic table names
  expect_true(file.exists(db_path))
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_true("reference_final" %in% tables)
  expect_true("qc_reference" %in% tables)
})


# ------------------------------------------------------------
# Ambiguity checking
# ------------------------------------------------------------

test_that("rb_curate_reference stops on unresolved ambiguity by default", {
  
  ambiguous_input <- data.frame(
    species = c("Species A", "Species B"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT"),
    seq_type = "COI",
    source = "test",
    stringsAsFactors = FALSE
  )
  
  tax_table <- make_wrapper_taxonomy(unique(ambiguous_input$species))
  
  expect_error(
    rb_curate_reference(
      ambiguous_input,
      taxonomy_table = tax_table,
      check_ambiguity = TRUE,
      on_ambiguous = "stop"
    ),
    "ambiguous"
  )
})


test_that("rb_curate_reference warns on unresolved ambiguity when on_ambiguous = 'warn'", {
  
  ambiguous_input <- data.frame(
    species = c("Species A", "Species B"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT"),
    seq_type = "COI",
    source = "test",
    stringsAsFactors = FALSE
  )
  
  tax_table <- make_wrapper_taxonomy(unique(ambiguous_input$species))
  
  expect_warning(
    rb_curate_reference(
      ambiguous_input,
      taxonomy_table = tax_table,
      check_ambiguity = TRUE,
      on_ambiguous = "warn",
      db_path = NULL
    ),
    "ambiguous"
  )
})


test_that("rb_curate_reference ignores ambiguity when on_ambiguous = 'ignore'", {
  
  ambiguous_input <- data.frame(
    species = c("Species A", "Species B"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT"),
    seq_type = "COI",
    source = "test",
    stringsAsFactors = FALSE
  )
  
  tax_table <- make_wrapper_taxonomy(unique(ambiguous_input$species))
  
  expect_no_error(
    rb_curate_reference(
      ambiguous_input,
      taxonomy_table = tax_table,
      check_ambiguity = TRUE,
      on_ambiguous = "ignore",
      db_path = NULL
    )
  )
})


test_that("rb_curate_reference skips ambiguity check when check_ambiguity = FALSE", {
  
  ambiguous_input <- data.frame(
    species = c("Species A", "Species B"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT"),
    seq_type = "COI",
    source = "test",
    stringsAsFactors = FALSE
  )
  
  tax_table <- make_wrapper_taxonomy(unique(ambiguous_input$species))
  
  expect_no_error(
    rb_curate_reference(
      ambiguous_input,
      taxonomy_table = tax_table,
      check_ambiguity = FALSE,
      db_path = NULL
    )
  )
})


# ------------------------------------------------------------
# Taxonomy validation
# ------------------------------------------------------------

test_that("rb_curate_reference errors when taxonomy is incomplete and require_taxonomy = TRUE", {
  
  input_data <- make_wrapper_test_data()
  
  # Taxonomy table missing genus
  incomplete_tax <- data.frame(
    species = unique(input_data$species),
    family = "Testidae",
    order = "Testorder",
    class = "Testclass",
    phylum = "Testphylum",
    kingdom = "Testkingdom",
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_curate_reference(
      data = input_data,
      taxonomy_table = incomplete_tax,
      require_taxonomy = TRUE,
      db_path = NULL
    ),
    "taxonomy"
  )
})


test_that("rb_curate_reference allows incomplete taxonomy when require_taxonomy = FALSE", {
  
  input_data <- make_wrapper_test_data()
  
  incomplete_tax <- data.frame(
    species = unique(input_data$species),
    family = "Testidae",
    stringsAsFactors = FALSE
  )
  
  expect_no_error(
    rb_curate_reference(
      data = input_data,
      taxonomy_table = incomplete_tax,
      require_taxonomy = FALSE,
      db_path = NULL
    )
  )
})


# ------------------------------------------------------------
# ML classification branch
# ------------------------------------------------------------

test_that("rb_curate_reference runs ML classification and updates seq_type", {
  
  skip_if_not_installed("parsnip")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("ranger")
  skip_on_cran()
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    run_classifier = TRUE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE,
    classifier_trees = 100,
    db_path = NULL
  )
  
  expect_type(result, "list")
  
  # Classifier should have trained
  expect_false(is.null(result$classifier))
  
  # In qc1, we had 5 "other" sequences
  expect_equal(sum(result$qc1$seq_type == "other"), 5)
  
  # In qc2, at least some "other" sequences should have been reclassified
  expect_true(sum(result$qc2$seq_type == "other") < 5)
  
  # Temporary row ID should be removed
  expect_false(".curate_row_id" %in% names(result$final_data))
})


# ------------------------------------------------------------
# Barcode gap analysis branch
# ------------------------------------------------------------

test_that("rb_curate_reference can run optional barcode gap analysis", {
  
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  skip_on_cran()
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = TRUE,
    barcode_gap_min_seqs = 2,
    barcode_gap_parallel = FALSE,
    db_path = NULL
  )
  
  expect_true("barcode_gap" %in% names(result))
})


test_that("rb_curate_reference skips barcode gap when run_barcode_gap = FALSE", {
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE,
    db_path = NULL
  )
  
  expect_true("barcode_gap" %in% names(result))
  expect_null(result$barcode_gap)
})


# ------------------------------------------------------------
# Ambiguity archiving
# ------------------------------------------------------------

test_that("rb_curate_reference archives ambiguity table if present", {
  
  skip_if_not_installed("RSQLite")
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  
  tmp_dir <- tempfile("barcurateR_ambiguity_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  ambiguity_path <- file.path(tmp_dir, "ambiguous_sequences.csv")
  
  ambiguity_table <- data.frame(
    species = "Species_1",
    sequence = "ACGTACGTACGT",
    note = "test ambiguity",
    stringsAsFactors = FALSE
  )
  
  utils::write.csv(ambiguity_table, ambiguity_path, row.names = FALSE)
  
  db_path <- file.path(tmp_dir, "test.db")
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    db_path = db_path,
    ambiguity_path = ambiguity_path,
    archive_ambiguity = TRUE,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE
  )
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_true("ambiguous_sequences" %in% tables)
})


test_that("rb_curate_reference skips ambiguity archiving when archive_ambiguity = FALSE", {
  
  skip_if_not_installed("RSQLite")
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  
  tmp_dir <- tempfile("barcurateR_no_ambiguity_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  ambiguity_path <- file.path(tmp_dir, "ambiguous_sequences.csv")
  
  ambiguity_table <- data.frame(
    species = "Species_1",
    sequence = "ACGTACGTACGT",
    stringsAsFactors = FALSE
  )
  
  utils::write.csv(ambiguity_table, ambiguity_path, row.names = FALSE)
  
  db_path <- file.path(tmp_dir, "test.db")
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    db_path = db_path,
    ambiguity_path = ambiguity_path,
    archive_ambiguity = FALSE,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE
  )
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_false("ambiguous_sequences" %in% tables)
})


# ------------------------------------------------------------
# Custom table names
# ------------------------------------------------------------

test_that("rb_curate_reference respects custom table names", {
  
  skip_if_not_installed("RSQLite")
  
  input_data <- make_wrapper_test_data()
  tax_table <- make_wrapper_taxonomy(unique(input_data$species))
  db_path <- tempfile(fileext = ".db")
  
  result <- rb_curate_reference(
    data = input_data,
    taxonomy_table = tax_table,
    db_path = db_path,
    reference_table = "my_reference",
    qc_table = "my_qc",
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE
  )
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_true("my_reference" %in% tables)
  expect_true("my_qc" %in% tables)
  expect_false("reference_final" %in% tables)
  expect_false("qc_reference" %in% tables)
})


# ------------------------------------------------------------
# Missing required columns
# ------------------------------------------------------------

test_that("rb_curate_reference handles missing required columns gracefully", {
  
  bad_data <- data.frame(
    species = "Species A",
    sequence = "ACGT",
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_curate_reference(bad_data),
    "Missing required columns"
  )
})