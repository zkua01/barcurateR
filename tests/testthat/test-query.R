# tests/testthat/test-query.R

test_that("rb_get_sequences retrieves and filters sequences correctly", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Default retrieval (qc_flag = "pass" by default)
  res <- rb_get_sequences(con)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5) # 1 contaminant excluded
  
  # Filter by marker
  res_coi <- rb_get_sequences(con, marker = "COI")
  expect_equal(nrow(res_coi), 2)
  
  # Filter by species
  res_sp <- rb_get_sequences(con, species = "Species A")
  expect_equal(nrow(res_sp), 2)
  
  # Include all qc_flags by setting to NULL
  res_all <- rb_get_sequences(con, qc_flag = NULL)
  expect_equal(nrow(res_all), 6)
})

test_that("rb_get_sequences uses custom table_name", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Write to a custom table
  DBI::dbWriteTable(con, "my_custom_ref", DBI::dbReadTable(con, "reference_final"))
  
  res <- rb_get_sequences(con, table_name = "my_custom_ref", qc_flag = NULL)
  expect_equal(nrow(res), 6)
})

test_that("rb_list_taxa lists taxa at different ranks", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res_species <- rb_list_taxa(con, rank = "species")
  expect_equal(nrow(res_species), 3)
  
  res_genus <- rb_list_taxa(con, rank = "genus")
  expect_equal(nrow(res_genus), 3)
  
  res_family <- rb_list_taxa(con, rank = "family")
  expect_equal(nrow(res_family), 2)
})

test_that("rb_filter_reference filters in-memory data frames", {
  df <- data.frame(
    species = c("A", "B", "C"),
    seq_type = c("COI", "12S", "COI"),
    qc_flag = c("pass", "pass", "contaminant"),
    sequence = c("ACGT", "TTTT", "GGGG"),
    stringsAsFactors = FALSE
  )
  
  res <- rb_filter_reference(df, marker = "COI", qc_flag = "pass")
  expect_equal(nrow(res), 1)
  expect_equal(res$species, "A")
})

test_that("rb_taxonomy retrieves taxonomy and handles optional columns", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  res <- rb_taxonomy(con)
  
  # Required ranks must be present
  expect_true(all(c("kingdom", "phylum", "class", "order", "family", "genus", "species") %in% names(res)))
  
  # Optional columns should be present because they exist in the mock DB
  expect_true("occurrence" %in% names(res))
  expect_true("habitat" %in% names(res))
  
  # Filter by species
  res_sp <- rb_taxonomy(con, species = "Species A")
  expect_equal(nrow(res_sp), 1)
})

test_that("rb_taxonomy works when optional columns are missing", {
  con <- create_mock_db()
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Create a table without occurrence/habitat
  ref_no_eco <- DBI::dbReadTable(con, "reference_final")
  ref_no_eco$occurrence <- NULL
  ref_no_eco$habitat <- NULL
  DBI::dbWriteTable(con, "ref_no_eco", ref_no_eco)
  
  res <- rb_taxonomy(con, table_name = "ref_no_eco")
  
  expect_true(all(c("kingdom", "species") %in% names(res)))
  expect_false("occurrence" %in% names(res))
  expect_false("habitat" %in% names(res))
})