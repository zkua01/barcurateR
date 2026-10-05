# tests/testthat/test-connect.R

# ============================================================
# Tests for connect.R functions: rb_connect, rb_disconnect, 
# rb_tables, rb_read_table, and rb_available_datasets.
# ============================================================

test_that("rb_available_datasets returns a properly formatted data frame", {
  skip_if_not_installed("withr")
  skip_if_not_installed("RSQLite")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-avail-test")
  unlink(cache_dir, recursive = TRUE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Create real SQLite DB files in cache
  fake_full_db <- file.path(cache_dir, "YZFishDB.db")
  con1 <- DBI::dbConnect(RSQLite::SQLite(), fake_full_db)
  DBI::dbWriteTable(con1, "test_table", data.frame(x = 1))
  DBI::dbDisconnect(con1)
  
  fake_user_db <- file.path(cache_dir, "my_custom_db.sqlite")
  con2 <- DBI::dbConnect(RSQLite::SQLite(), fake_user_db)
  DBI::dbWriteTable(con2, "test_table", data.frame(x = 1))
  DBI::dbDisconnect(con2)
  
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir
  )
  
  datasets <- rb_available_datasets()
  
  expect_true(is.data.frame(datasets))
  expect_true(all(c("name", "version", "path") %in% names(datasets)))
  
  # Should find our fake databases
  expect_true(any(grepl("YZFishDB", datasets$path)))
  expect_true(any(grepl("my_custom_db", datasets$path)))
  
  # Check versions
  yz_row <- datasets[grepl("YZFishDB", datasets$path), ]
  expect_equal(yz_row$version, "full-local")
  
  user_row <- datasets[grepl("my_custom_db", datasets$path), ]
  expect_equal(user_row$version, "user-curated")
})

test_that("rb_connect prioritizes full-local YZFishDB over other databases", {
  skip_if_not_installed("withr")
  skip_if_not_installed("RSQLite")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-connect-priority-test")
  unlink(cache_dir, recursive = TRUE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Create a valid SQLite DB as user-curated
  user_db <- file.path(cache_dir, "user_db.sqlite")
  con_user <- DBI::dbConnect(RSQLite::SQLite(), user_db)
  DBI::dbWriteTable(con_user, "test_table", data.frame(x = 1))
  DBI::dbDisconnect(con_user)
  
  # Create a valid SQLite DB as full-local (YZFishDB.db)
  full_db <- file.path(cache_dir, "YZFishDB.db")
  con_full <- DBI::dbConnect(RSQLite::SQLite(), full_db)
  DBI::dbWriteTable(con_full, "yzfishdb_table", data.frame(x = 2))
  DBI::dbDisconnect(con_full)
  
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir
  )
  
  # Connect without path
  con <- rb_connect()
  on.exit(rb_disconnect(con), add = TRUE)
  
  expect_true(DBI::dbIsValid(con))
  
  # Check which database was opened by looking at its tables
  tables <- DBI::dbListTables(con)
  expect_true("yzfishdb_table" %in% tables)
  expect_false("test_table" %in% tables)
})

test_that("rb_connect falls back to user-curated if no full-local exists", {
  skip_if_not_installed("withr")
  skip_if_not_installed("RSQLite")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-connect-fallback-test")
  unlink(cache_dir, recursive = TRUE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Only create a user-curated DB
  user_db <- file.path(cache_dir, "my_legacy_regionbarcoder.sqlite")
  con_user <- DBI::dbConnect(RSQLite::SQLite(), user_db)
  DBI::dbWriteTable(con_user, "legacy_table", data.frame(x = 1))
  DBI::dbDisconnect(con_user)
  
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir
  )
  
  # Mock rb_available_datasets to ONLY return the user DB so we don't accidentally 
  # find the actual package's bundled extdata database (which might be YZFishDB or demo)
  testthat::local_mocked_bindings(
    rb_available_datasets = function() {
      data.frame(
        name = "my_legacy_regionbarcoder",
        version = "user-curated",
        path = user_db,
        stringsAsFactors = FALSE
      )
    },
    .package = "barcurateR"
  )
  
  con <- rb_connect()
  on.exit(rb_disconnect(con), add = TRUE)
  
  expect_true(DBI::dbIsValid(con))
  tables <- DBI::dbListTables(con)
  expect_true("legacy_table" %in% tables)
})

test_that("rb_connect errors when no local databases are found and download is FALSE", {
  testthat::local_mocked_bindings(
    rb_available_datasets = function() {
      data.frame(name = character(0), version = character(0), path = character(0), stringsAsFactors = FALSE)
    },
    .package = "barcurateR"
  )
  
  expect_error(rb_connect(download = FALSE), "No local databases found")
})

test_that("rb_tables and rb_read_table work correctly", {
  skip_if_not_installed("RSQLite")
  
  con <- create_mock_db() # Uses the helper-mock-db.R function
  on.exit(rb_disconnect(con), add = TRUE)
  
  tables <- rb_tables(con)
  expect_true(is.data.frame(tables))
  expect_true(all(c("table", "n_rows") %in% names(tables)))
  expect_true("reference_final" %in% tables$table)
  expect_equal(tables$n_rows[tables$table == "reference_final"], 6)
  
  refs <- rb_read_table(con, "reference_final")
  expect_equal(nrow(refs), 6)
  expect_true("species" %in% names(refs))
  
  expect_error(rb_read_table(con, "non_existent_table"), "Table not found")
})

test_that("rb_disconnect handles invalid connections gracefully", {
  con <- create_mock_db()
  DBI::dbDisconnect(con)
  
  # Should not error even if already disconnected
  expect_invisible(rb_disconnect(con))
  expect_true(rb_disconnect(con))
})