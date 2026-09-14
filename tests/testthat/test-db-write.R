# ============================================================
# tests/testthat/test-db-write.R
#
# Tests for SQLite writing functions.
#
# Updated to use generic table names:
#   reference_final
#   qc_reference
#   barcode_gap_metrics
#   ambiguous_sequences
#
# Taxonomy (kingdom-to-species) is required.
# occurrence and habitat are optional.
# ============================================================


# ------------------------------------------------------------
# Helper: test data with full taxonomy
# ------------------------------------------------------------

make_reference_data <- function() {
  data.frame(
    species = c("Danio rerio", "Danio rerio", "Cyprinus carpio"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT", "TTTTGGGGCCCC"),
    seq_type = c("COI", "COI", "12S"),
    source = c("ncbi", "ncbi", "bold"),
    qc_flag = c("pass", "pass", "pass"),
    unique_code = c("seq_00001", "seq_00001", "seq_00002"),
    kingdom = "Animalia",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Cypriniformes",
    family = "Cyprinidae",
    genus = c("Danio", "Danio", "Cyprinus"),
    stringsAsFactors = FALSE
  )
}

make_qc_data <- function() {
  data.frame(
    unique_code = c("seq_00001", "seq_00002"),
    species = c("Danio rerio", "Cyprinus carpio"),
    qc_flag = c("pass", "contaminant"),
    stringsAsFactors = FALSE
  )
}

make_barcode_gap_data <- function() {
  data.frame(
    species = c("Danio rerio", "Cyprinus carpio"),
    marker = c("COI", "12S"),
    gap_exists = c(TRUE, FALSE),
    gap_width = c(0.05, NA_real_),
    stringsAsFactors = FALSE
  )
}

make_ambiguity_data <- function() {
  data.frame(
    species = c("Danio rerio", "Carassius auratus"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT"),
    note = c("same sequence, two species", "same sequence, two species"),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# rb_prepare_for_sqlite()
# ------------------------------------------------------------

test_that("rb_prepare_for_sqlite converts factors to character", {
  
  df <- data.frame(
    species = factor("Danio rerio"),
    source = factor("ncbi"),
    length = 100,
    stringsAsFactors = FALSE
  )
  
  out <- rb_prepare_for_sqlite(df)
  
  expect_type(out$species, "character")
  expect_type(out$source, "character")
})


test_that("rb_prepare_for_sqlite rejects list-columns", {
  
  df <- data.frame(
    species = "Danio rerio",
    sequence = "ACGT",
    stringsAsFactors = FALSE
  )
  
  df$bad_column <- list(1:3)
  
  expect_error(
    rb_prepare_for_sqlite(df),
    "list column"
  )
})


# ------------------------------------------------------------
# rb_prepare_reference_table()
# ------------------------------------------------------------

test_that("rb_prepare_reference_table ensures required columns", {
  
  df <- data.frame(
    species = "Danio rerio",
    sequence = "ACGTACGTACGT",
    kingdom = "Animalia",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Cypriniformes",
    family = "Cyprinidae",
    genus = "Danio",
    stringsAsFactors = FALSE
  )
  
  out <- rb_prepare_reference_table(df, require_taxonomy = TRUE)
  
  expect_true("qc_flag" %in% names(out))
  expect_true("source" %in% names(out))
  expect_true("seq_type" %in% names(out))
  expect_true("unique_code" %in% names(out))
})


test_that("rb_prepare_reference_table requires full taxonomy when require_taxonomy = TRUE", {
  
  df <- data.frame(
    species = "Danio rerio",
    sequence = "ACGTACGTACGT",
    kingdom = "Animalia",
    phylum = "Chordata",
    class = "Actinopterygii",
    order = "Cypriniformes",
    family = "Cyprinidae",
    # genus is missing
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_prepare_reference_table(df, require_taxonomy = TRUE),
    "taxonomy"
  )
})


test_that("rb_prepare_reference_table skips taxonomy check when require_taxonomy = FALSE", {
  
  df <- data.frame(
    species = "Danio rerio",
    sequence = "ACGTACGTACGT",
    stringsAsFactors = FALSE
  )
  
  expect_no_error(
    rb_prepare_reference_table(df, require_taxonomy = FALSE)
  )
})


test_that("occurrence and habitat are optional", {
  
  df_no_ecology <- make_reference_data()
  
  expect_no_error(
    rb_prepare_reference_table(df_no_ecology, require_taxonomy = TRUE)
  )
  
  df_with_ecology <- make_reference_data()
  df_with_ecology$occurrence <- "native"
  df_with_ecology$habitat <- "freshwater"
  
  expect_no_error(
    rb_prepare_reference_table(df_with_ecology, require_taxonomy = TRUE)
  )
})


# ------------------------------------------------------------
# rb_write_reference_table()
# ------------------------------------------------------------

test_that("rb_write_reference_table writes to reference_final by default", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  data <- make_reference_data()
  
  rb_write_reference_table(con, data)
  
  expect_true("reference_final" %in% DBI::dbListTables(con))
  
  out <- DBI::dbReadTable(con, "reference_final")
  expect_equal(nrow(out), 3)
})


test_that("rb_write_reference_table can write to a custom table name", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  data <- make_reference_data()
  
  rb_write_reference_table(con, data, table_name = "custom_reference")
  
  expect_true("custom_reference" %in% DBI::dbListTables(con))
  expect_false("reference_final" %in% DBI::dbListTables(con))
})


test_that("rb_write_reference_table errors without full taxonomy", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  bad_data <- make_reference_data()
  bad_data$genus <- NULL
  
  expect_error(
    rb_write_reference_table(con, bad_data, require_taxonomy = TRUE),
    "taxonomy"
  )
})


# ------------------------------------------------------------
# rb_write_qc_table()
# ------------------------------------------------------------

test_that("rb_write_qc_table writes to qc_reference by default", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  qc_data <- make_qc_data()
  
  rb_write_qc_table(con, qc_data)
  
  expect_true("qc_reference" %in% DBI::dbListTables(con))
  
  out <- DBI::dbReadTable(con, "qc_reference")
  expect_equal(nrow(out), 2)
})


test_that("rb_write_qc_table adds qc_flag if missing", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  qc_no_flag <- data.frame(
    unique_code = "seq_00001",
    species = "Danio rerio",
    stringsAsFactors = FALSE
  )
  
  rb_write_qc_table(con, qc_no_flag)
  
  out <- DBI::dbReadTable(con, "qc_reference")
  expect_true("qc_flag" %in% names(out))
  expect_equal(out$qc_flag, "pass")
})


# ------------------------------------------------------------
# rb_write_barcode_gap_table()
# ------------------------------------------------------------

test_that("rb_write_barcode_gap_table writes to barcode_gap_metrics by default", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  gap_data <- make_barcode_gap_data()
  
  rb_write_barcode_gap_table(con, gap_data)
  
  expect_true("barcode_gap_metrics" %in% DBI::dbListTables(con))
  
  out <- DBI::dbReadTable(con, "barcode_gap_metrics")
  expect_equal(nrow(out), 2)
})


test_that("rb_write_barcode_gap_table requires species and marker columns", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  bad_gap <- data.frame(x = 1)
  
  expect_error(
    rb_write_barcode_gap_table(con, bad_gap),
    "Missing required columns"
  )
})


# ------------------------------------------------------------
# rb_write_ambiguity_table()
# ------------------------------------------------------------

test_that("rb_write_ambiguity_table writes to ambiguous_sequences by default", {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  amb_data <- make_ambiguity_data()
  
  rb_write_ambiguity_table(con, amb_data)
  
  expect_true("ambiguous_sequences" %in% DBI::dbListTables(con))
  
  out <- DBI::dbReadTable(con, "ambiguous_sequences")
  expect_equal(nrow(out), 2)
})


# ------------------------------------------------------------
# rb_build_reference_db()
# ------------------------------------------------------------

test_that("rb_build_reference_db creates a SQLite file with all tables", {
  
  db_path <- tempfile(fileext = ".db")
  
  rb_build_reference_db(
    path = db_path,
    reference_data = make_reference_data(),
    qc_data = make_qc_data(),
    barcode_gap_data = make_barcode_gap_data(),
    ambiguity_data = make_ambiguity_data()
  )
  
  expect_true(file.exists(db_path))
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_true("reference_final" %in% tables)
  expect_true("qc_reference" %in% tables)
  expect_true("barcode_gap_metrics" %in% tables)
  expect_true("ambiguous_sequences" %in% tables)
})


test_that("rb_build_reference_db works with only reference data", {
  
  db_path <- tempfile(fileext = ".db")
  
  rb_build_reference_db(
    path = db_path,
    reference_data = make_reference_data()
  )
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_true("reference_final" %in% tables)
  expect_false("qc_reference" %in% tables)
  expect_false("barcode_gap_metrics" %in% tables)
  expect_false("ambiguous_sequences" %in% tables)
})


test_that("rb_build_reference_db respects custom table names", {
  
  db_path <- tempfile(fileext = ".db")
  
  rb_build_reference_db(
    path = db_path,
    reference_data = make_reference_data(),
    qc_data = make_qc_data(),
    reference_table = "my_reference",
    qc_table = "my_qc"
  )
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  tables <- DBI::dbListTables(con)
  DBI::dbDisconnect(con)
  
  expect_true("my_reference" %in% tables)
  expect_true("my_qc" %in% tables)
  expect_false("reference_final" %in% tables)
  expect_false("qc_reference" %in% tables)
})


test_that("rb_build_reference_db errors without taxonomy when require_taxonomy = TRUE", {
  
  db_path <- tempfile(fileext = ".db")
  
  bad_data <- make_reference_data()
  bad_data$genus <- NULL
  
  expect_error(
    rb_build_reference_db(
      path = db_path,
      reference_data = bad_data,
      require_taxonomy = TRUE
    ),
    "taxonomy"
  )
})


test_that("rb_build_reference_db allows missing taxonomy when require_taxonomy = FALSE", {
  
  db_path <- tempfile(fileext = ".db")
  
  minimal_data <- data.frame(
    species = "Danio rerio",
    sequence = "ACGTACGTACGT",
    stringsAsFactors = FALSE
  )
  
  expect_no_error(
    rb_build_reference_db(
      path = db_path,
      reference_data = minimal_data,
      require_taxonomy = FALSE
    )
  )
})