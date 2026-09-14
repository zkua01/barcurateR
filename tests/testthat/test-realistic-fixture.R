# ============================================================
# tests/testthat/test-realistic-fixture.R
#
# End-to-end integration test using realistic fixture data:
#
#   - NCBI marker sequences
#   - BOLD COI sequences
#   - Quick taxonomy table
#
# This test assumes the user has already prepared the input
# tables in a format compatible with rb_parse_source_table().
# ============================================================


# ------------------------------------------------------------
# Helper: clean NCBI fixture
#
# The downloaded NCBI fixture contains some records that are
# not ideal reference sequences, for example:
#
#   - PREDICTED mRNA records
#   - duplicate accessions assigned to multiple markers
#
# This helper performs light user-side cleanup before
# standardization.
# ------------------------------------------------------------

clean_ncbi_fixture <- function(raw) {
  
  # Remove predicted / non-reference records
  bad_pattern <- paste0(
    "PREDICTED|mRNA|synthetic construct|",
    "uncultured|environmental sample|vector"
  )
  
  raw <- raw[
    !grepl(bad_pattern, raw$description, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  
  # Keep only one record per accession.
  # This avoids the same multi-gene accession being treated as
  # both 12S and 16S.
  raw <- raw[!duplicated(raw$sequence_id), , drop = FALSE]
  
  # Use the marker query as the marker description for parsing.
  # This gives rb_extract_markers() a clean marker label such as
  # "12S", "16S", or "COI".
  raw$description <- raw$marker_query
  
  raw
}


# ------------------------------------------------------------
# Main realistic fixture test
# ------------------------------------------------------------

test_that("realistic fixture can be standardized, curated, and written to SQLite", {
  
  skip_if_not_installed("RSQLite")
  skip_on_cran()
  
  # ----------------------------------------------------------
  # 1. Locate fixture files
  # ----------------------------------------------------------
  
  ncbi_path <- testthat::test_path("fixtures", "ncbi_marker_sequences.csv")
  bold_path <- testthat::test_path("fixtures", "bold_coi_sequences.tsv")
  tax_path <- testthat::test_path("fixtures", "quick_taxonomy.csv")
  
  skip_if_not(file.exists(ncbi_path), "NCBI fixture missing")
  skip_if_not(file.exists(bold_path), "BOLD fixture missing")
  skip_if_not(file.exists(tax_path), "Taxonomy fixture missing")
  
  # ----------------------------------------------------------
  # 2. Load fixture files
  # ----------------------------------------------------------
  
  ncbi_raw <- readr::read_csv(
    ncbi_path,
    show_col_types = FALSE
  )
  
  bold_raw <- readr::read_tsv(
    bold_path,
    show_col_types = FALSE
  )
  
  taxonomy_table <- readr::read_csv(
    tax_path,
    show_col_types = FALSE
  )
  
  expect_true(nrow(ncbi_raw) > 0)
  expect_true(nrow(bold_raw) > 0)
  expect_true(nrow(taxonomy_table) > 0)
  
  # ----------------------------------------------------------
  # 3. Light user-side cleanup of NCBI fixture
  # ----------------------------------------------------------
  
  ncbi_raw <- clean_ncbi_fixture(ncbi_raw)
  
  expect_true(nrow(ncbi_raw) > 0)
  
  # ----------------------------------------------------------
  # 4. Prepare BOLD fixture for standardization
  #
  # The current standardization function expects a description
  # column from which the marker can be extracted. For BOLD,
  # we can use the marker label directly.
  # ----------------------------------------------------------
  
  bold_raw$description <- "COI"
  
  # ----------------------------------------------------------
  # 5. Standardize sources
  # ----------------------------------------------------------
  
  ncbi_parsed <- rb_parse_source_table(
    ncbi_raw,
    source_name = "ncbi",
    column_map = c(
      sequence_id = "sequence_id",
      species = "species_query",
      sequence = "sequence"
    ),
    description_col = "description"
  )
  
  bold_parsed <- rb_parse_source_table(
    bold_raw,
    source_name = "bold",
    column_map = c(
      sequence_id = "processid",
      species = "species",
      sequence = "nuc"
    ),
    description_col = "description"
  )
  
  expect_true(nrow(ncbi_parsed) > 0)
  expect_true(nrow(bold_parsed) > 0)
  
  expect_true(all(c("species", "sequence", "seq_type") %in% names(ncbi_parsed)))
  expect_true(all(c("species", "sequence", "seq_type") %in% names(bold_parsed)))
  
  # The fixture should contain at least COI and one rRNA marker
  expect_true(any(ncbi_parsed$seq_type == "COI"))
  expect_true(
    any(ncbi_parsed$seq_type %in% c("12S", "16S"))
  )
  
  # BOLD fixture should be COI
  expect_true(all(bold_parsed$seq_type == "COI"))
  
  # ----------------------------------------------------------
  # 6. Combine sources
  # ----------------------------------------------------------
  
  combined <- rb_combine_sources(
    list(
      ncbi_parsed,
      bold_parsed
    ),
    species_col = "species",
    sequence_col = "sequence"
  )
  
  expect_true(nrow(combined) > 0)
  
  # ----------------------------------------------------------
  # 7. Resolve ambiguities before QC
  # ----------------------------------------------------------
  
  combined <- rb_resolve_ambiguous(
    combined,
    on_unresolved = "drop"
  )
  
  expect_true(nrow(combined) > 0)
  
  # ----------------------------------------------------------
  # 8. Run the high-level curation wrapper
  # ----------------------------------------------------------
  
  db_path <- tempfile(fileext = ".db")
  
  result <- rb_curate_reference(
    data = combined,
    taxonomy_table = taxonomy_table,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE,
    db_path = db_path,
    check_ambiguity = TRUE,
    on_ambiguous = "stop",
    require_taxonomy = TRUE
  )
  
  # ----------------------------------------------------------
  # 9. Validate wrapper output
  # ----------------------------------------------------------
  
  expect_type(result, "list")
  
  expect_true(all(
    c("qc1", "qc2", "final_data", "db_path") %in% names(result)
  ))
  
  expect_true(nrow(result$final_data) > 0)
  
  # QC flags should be compiled
  expect_true("qc_flag" %in% names(result$qc1))
  expect_true("qc_flag" %in% names(result$qc2))
  expect_true("qc_flag" %in% names(result$final_data))
  
  # Taxonomy should be joined
  expect_true("kingdom" %in% names(result$final_data))
  expect_true("phylum" %in% names(result$final_data))
  expect_true("class" %in% names(result$final_data))
  expect_true("order" %in% names(result$final_data))
  expect_true("family" %in% names(result$final_data))
  expect_true("genus" %in% names(result$final_data))
  expect_true("species" %in% names(result$final_data))
  
  # Temporary row identifier should not be exported
  expect_false(".curate_row_id" %in% names(result$final_data))
  expect_false(".curate_row_id" %in% names(result$qc2))
  
  # ----------------------------------------------------------
  # 10. Validate SQLite database
  # ----------------------------------------------------------
  
  expect_true(file.exists(db_path))
  
  con <- rb_connect(db_path)
  
  tables <- rb_tables(con)
  
  expect_s3_class(tables, "data.frame")
  expect_true("reference_final" %in% tables$table)
  expect_true("qc_reference" %in% tables$table)
  
  # ----------------------------------------------------------
  # 11. Legacy query compatibility
  # ----------------------------------------------------------
  
  pass_seqs <- rb_get_sequences(
    con,
    qc_flag = "pass"
  )
  
  expect_s3_class(pass_seqs, "data.frame")
  
  marker_cov <- rb_marker_coverage(con)
  
  expect_s3_class(marker_cov, "data.frame")
  
  qc_sum <- rb_qc_summary(con)
  
  expect_s3_class(qc_sum, "data.frame")
  expect_true(all(c("qc_flag", "n_sequences") %in% names(qc_sum)))
  
  tax <- rb_taxonomy(con)
  
  expect_s3_class(tax, "data.frame")
  
  expect_true(all(
    c(
      "kingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    ) %in% names(tax)
  ))
  
  rb_disconnect(con)
})


# ------------------------------------------------------------
# Optional export compatibility test
# ------------------------------------------------------------

test_that("realistic fixture output can be exported to DADA2 format", {
  
  skip_if_not_installed("RSQLite")
  skip_on_cran()
  
  ncbi_path <- testthat::test_path("fixtures", "ncbi_marker_sequences.csv")
  bold_path <- testthat::test_path("fixtures", "bold_coi_sequences.tsv")
  tax_path <- testthat::test_path("fixtures", "quick_taxonomy.csv")
  
  skip_if_not(file.exists(ncbi_path), "NCBI fixture missing")
  skip_if_not(file.exists(bold_path), "BOLD fixture missing")
  skip_if_not(file.exists(tax_path), "Taxonomy fixture missing")
  
  ncbi_raw <- readr::read_csv(
    ncbi_path,
    show_col_types = FALSE
  )
  
  bold_raw <- readr::read_tsv(
    bold_path,
    show_col_types = FALSE
  )
  
  taxonomy_table <- readr::read_csv(
    tax_path,
    show_col_types = FALSE
  )
  
  ncbi_raw <- clean_ncbi_fixture(ncbi_raw)
  bold_raw$description <- "COI"
  
  ncbi_parsed <- rb_parse_source_table(
    ncbi_raw,
    source_name = "ncbi",
    column_map = c(
      sequence_id = "sequence_id",
      species = "species_query",
      sequence = "sequence"
    ),
    description_col = "description"
  )
  
  bold_parsed <- rb_parse_source_table(
    bold_raw,
    source_name = "bold",
    column_map = c(
      sequence_id = "processid",
      species = "species",
      sequence = "nuc"
    ),
    description_col = "description"
  )
  
  combined <- rb_combine_sources(
    list(
      ncbi_parsed,
      bold_parsed
    ),
    species_col = "species",
    sequence_col = "sequence"
  )
  
  combined <- rb_resolve_ambiguous(
    combined,
    on_unresolved = "drop"
  )
  
  result <- rb_curate_reference(
    data = combined,
    taxonomy_table = taxonomy_table,
    run_classifier = FALSE,
    run_divergence = FALSE,
    run_barcode_gap = FALSE,
    db_path = NULL,
    check_ambiguity = TRUE,
    on_ambiguous = "stop",
    require_taxonomy = TRUE
  )
  
  # Use passing sequences if available; otherwise use final data
  export_data <- result$final_data[
    result$final_data$qc_flag == "pass",
    ,
    drop = FALSE
  ]
  
  if (nrow(export_data) == 0) {
    export_data <- result$final_data
  }
  
  expect_true(nrow(export_data) > 0)
  
  out_fasta <- tempfile(fileext = ".fasta")
  
  rb_export_dada2(
    export_data,
    out_fasta
  )
  
  expect_true(file.exists(out_fasta))
  
  # Basic FASTA sanity check
  fasta_lines <- readLines(out_fasta, warn = FALSE)
  
  expect_true(any(grepl("^>", fasta_lines)))
})