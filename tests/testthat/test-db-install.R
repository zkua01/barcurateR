# ============================================================
# tests/testthat/test-db-install.R
#
# Tests for database installation, caching, and on-demand
# download/install behavior.
#
# These tests avoid network access by using local files as
# "download" sources.
# ============================================================


test_that("rb_db_dir uses environment override", {
  
  skip_if_not_installed("withr")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-env-test")
  unlink(cache_dir, recursive = TRUE)
  
  # Set both old and new environment names for compatibility
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir,
    BC_YZFISHDB_CACHE_DIR = cache_dir
  )
  
  observed <- rb_db_dir()
  
  expect_identical(
    observed,
    normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
  )
})


test_that("rb_install_db downloads to the cache and reuses existing files", {
  
  skip_if_not_installed("withr")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-download-test")
  unlink(cache_dir, recursive = TRUE)
  
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir,
    BC_YZFISHDB_CACHE_DIR = cache_dir
  )
  
  source_db <- file.path(tempdir(), "source-payload.db")
  writeLines("fake sqlite payload", source_db)
  
  installed <- rb_install_db(
    url = source_db,
    quiet = TRUE
  )
  
  expect_true(file.exists(installed))
  expect_equal(readLines(installed), "fake sqlite payload")
  
  # Change the source payload. If caching works, the installed file
  # should not be overwritten because overwrite = FALSE by default.
  writeLines("replacement payload", source_db)
  
  reused <- rb_install_db(
    url = source_db,
    quiet = TRUE
  )
  
  expect_identical(reused, installed)
  expect_equal(readLines(reused), "fake sqlite payload")
})


test_that("rb_install_db temporarily increases timeout for large downloads", {
  
  skip_if_not_installed("withr")
  
  observed <- new.env()
  observed$timeout <- NULL
  
  fake_downloader <- function(url, destfile, quiet = FALSE) {
    observed$timeout <- getOption("timeout")
    writeLines("timeout payload", destfile)
    invisible(destfile)
  }
  
  source_db <- file.path(tempdir(), "timeout-source.db")
  target_db <- file.path(tempdir(), "timeout-target.db")
  
  writeLines("timeout source", source_db)
  unlink(target_db)
  
  withr::local_options(
    barcurateR.download_file = fake_downloader
  )
  
  old_timeout <- getOption("timeout")
  
  rb_install_db(
    url = source_db,
    destfile = target_db,
    timeout = 3600,
    quiet = TRUE
  )
  
  expected_timeout <- max(
    3600,
    as.numeric(old_timeout),
    na.rm = TRUE
  )
  
  expect_equal(observed$timeout, expected_timeout)
  
  # Timeout option should be restored after installation
  expect_equal(getOption("timeout"), old_timeout)
})


test_that("rb_install_db checks md5 when requested", {
  
  source_db <- file.path(tempdir(), "md5-source.db")
  target_db <- file.path(tempdir(), "md5-target.db")
  
  writeLines("md5 payload", source_db)
  unlink(target_db)
  
  expect_error(
    rb_install_db(
      url = source_db,
      destfile = target_db,
      md5 = "not-the-right-md5",
      quiet = TRUE
    ),
    "md5"
  )
  
  # The target file should not be written if checksum validation fails
  expect_false(file.exists(target_db))
})


test_that("rb_install_db explains missing release URL when defaults are disabled", {
  
  skip_if_not_installed("withr")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-nourl-test")
  unlink(cache_dir, recursive = TRUE)
  
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir,
    BC_YZFISHDB_CACHE_DIR = cache_dir,
    BARCURATER_DB_URL = "",
    BC_YZFISHDB_URL = ""
  )
  
  withr::local_options(
    barcurateR.db_url = "",
    barcurateR.yzfishdb_url = ""
  )
  
  expect_error(
    rb_install_db(quiet = TRUE),
    "URL"
  )
})


test_that("rb_connect can install the database on demand", {
  
  skip_if_not_installed("withr")
  skip_if_not_installed("RSQLite")
  
  cache_dir <- file.path(tempdir(), "barcurateR-cache-connect-test")
  unlink(cache_dir, recursive = TRUE)
  
  # Create a real SQLite source database so rb_connect() can open it
  source_db <- file.path(tempdir(), "connect-source.sqlite")
  unlink(source_db)
  
  con_src <- DBI::dbConnect(RSQLite::SQLite(), source_db)
  DBI::dbWriteTable(
    con_src,
    "reference_final",
    data.frame(x = 1)
  )
  DBI::dbDisconnect(con_src)
  
  withr::local_envvar(
    BARCURATER_CACHE_DIR = cache_dir,
    BC_YZFISHDB_CACHE_DIR = cache_dir,
    BARCURATER_DB_URL = "",
    BC_YZFISHDB_URL = ""
  )
  
  withr::local_options(
    barcurateR.db_url = source_db,
    barcurateR.yzfishdb_url = source_db
  )
  
  con <- rb_connect(download = TRUE)
  
  expect_true(DBI::dbIsValid(con))
  expect_true("reference_final" %in% DBI::dbListTables(con))
  
  rb_disconnect(con)
})