#' Connect to and manage barcurateR SQLite databases
#'
#' @description
#' Functions to connect to, inspect, and manage barcurateR SQLite reference
#' databases. These functions provide the primary interface for interacting
#' with curated reference sequence databases.
#'
#' * `rb_connect()`: Opens a connection to a SQLite database.
#' * `rb_disconnect()`: Safely closes an active database connection.
#' * `rb_tables()`: Lists all tables in the database with row counts.
#' * `rb_read_table()`: Reads an entire table into a data frame.
#' * `rb_available_datasets()`: Lists all available bundled, cached, or legacy databases.
#'
#' @param path Path to a SQLite database file. If `NULL`, the default
#'   configured database is resolved (downloaded if requested, or the
#'   first available local database is used).
#' @param download Logical. If `TRUE` and the database file does not exist,
#'   attempt to download it from the default URL.
#' @param url The URL to download the database from (used when `download = TRUE`).
#' @param con A `DBIConnection` object returned by `rb_connect()`.
#' @param table The name of the table to read.
#'
#' @return
#' * `rb_connect()` returns a `DBIConnection` object (S4 class `SQLiteConnection`).
#' * `rb_disconnect()` returns `TRUE` invisibly.
#' * `rb_tables()` returns a data frame with columns `table` and `n_rows`.
#' * `rb_read_table()` returns a data frame containing the table contents.
#' * `rb_available_datasets()` returns a data frame with columns `name`, `version`, and `path`.
#'   It searches the package's `extdata` directory, the `barcurateR` cache directory,
#'   and the legacy `regionbarcoder` cache directory for any `.db` or `.sqlite` files.
#'
#' @details
#' The default database resolution order in `rb_connect()` when `path = NULL` is:
#' 1. If `download = TRUE` and no cached file exists, download from `rb_db_url()`.
#' 2. Otherwise, use the best available local database (preferring `full-local` YZFishDB
#'    over `demo` or `user-curated` databases).
#'
#' @examples
#' \dontrun{
#' # List all available local databases
#' avail <- rb_available_datasets()
#'
#' # Connect to the default database
#' con <- rb_connect()
#' rb_tables(con)
#'
#' # Read a specific table
#' refs <- rb_read_table(con, "reference_final")
#'
#' # Close the connection
#' rb_disconnect(con)
#' }
#'
#' @name rb_connect
#' @family database management
NULL

#' @rdname rb_connect
#' @export
rb_available_datasets <- function() {
  paths <- character(0)
  
  # 1. Check package extdata
  extdata_dir <- system.file("extdata", package = "barcurateR")
  if (nzchar(extdata_dir) && dir.exists(extdata_dir)) {
    ext_files <- list.files(extdata_dir, pattern = "\\.(sqlite|db)$", full.names = TRUE, ignore.case = TRUE)
    paths <- c(paths, ext_files)
  }
  
  # 2. Check barcurateR cache directory
  cache_dir <- tryCatch(rb_db_dir(create = FALSE), error = function(e) "")
  if (nzchar(cache_dir) && dir.exists(cache_dir)) {
    cache_files <- list.files(cache_dir, pattern = "\\.(sqlite|db)$", full.names = TRUE, ignore.case = TRUE)
    paths <- c(paths, cache_files)
  }
  
  # 3. Check legacy regionbarcoder cache directory
  legacy_dir <- tryCatch(tools::R_user_dir("regionbarcoder", which = "data"), error = function(e) "")
  if (nzchar(legacy_dir) && dir.exists(legacy_dir)) {
    legacy_files <- list.files(legacy_dir, pattern = "\\.(sqlite|db)$", full.names = TRUE, ignore.case = TRUE)
    paths <- c(paths, legacy_files)
  }
  
  paths <- unique(paths)
  
  if (length(paths) == 0) {
    return(data.frame(
      name = character(0),
      version = character(0),
      path = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  names_vec <- tools::file_path_sans_ext(basename(paths))
  versions <- ifelse(grepl("YZFishDB", basename(paths), ignore.case = TRUE), 
                     "full-local", 
                     ifelse(grepl("small_refdb|demo", basename(paths), ignore.case = TRUE),
                            "demo", 
                            "user-curated"))
  
  data.frame(
    name = names_vec,
    version = versions,
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

#' @rdname rb_connect
#' @export
rb_connect <- function(path = NULL, download = FALSE, url = rb_db_url()) {
  if (is.null(path)) {
    if (isTRUE(download) && !file.exists(rb_db_path())) {
      path <- rb_install_db(url = url)
    } else {
      avail <- rb_available_datasets()
      if (nrow(avail) == 0) {
        stop("No local databases found. Try setting `download = TRUE` or providing a `path`.", call. = FALSE)
      }
      # Prefer full-local YZFishDB, then fall back to demo or user-curated
      full_local_idx <- which(avail$version == "full-local")
      if (length(full_local_idx) > 0) {
        path <- avail$path[full_local_idx[1]]
      } else {
        path <- avail$path[1]
      }
    }
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Database file does not exist: ", path, call. = FALSE)
  }
  DBI::dbConnect(RSQLite::SQLite(), normalizePath(path, winslash = "/", mustWork = TRUE))
}

#' @rdname rb_connect
#' @export
rb_disconnect <- function(con) {
  if (DBI::dbIsValid(con)) {
    DBI::dbDisconnect(con)
  }
  invisible(TRUE)
}

#' @rdname rb_connect
#' @export
rb_tables <- function(con) {
  tables <- DBI::dbListTables(con)
  counts <- vapply(tables, function(tab) {
    quoted <- as.character(DBI::dbQuoteIdentifier(con, tab))
    DBI::dbGetQuery(con, paste0("select count(*) as n from ", quoted))[[1]]
  }, numeric(1))
  data.frame(table = tables, n_rows = as.integer(counts), stringsAsFactors = FALSE)
}

#' @rdname rb_connect
#' @export
rb_read_table <- function(con, table) {
  available <- DBI::dbListTables(con)
  if (!table %in% available) {
    stop("Table not found: ", table, call. = FALSE)
  }
  DBI::dbReadTable(con, table)
}

#' Internal helper to check if a file is a valid SQLite database
#' @noRd
rb_is_sqlite_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 16) {
    return(FALSE)
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  header <- readBin(con, "raw", n = 16)
  length(header) >= 16 && identical(header[1:15], charToRaw("SQLite format 3"))
}

#' @rdname rb_connect
#' @export
rb_available_datasets <- function() {
  paths <- character(0)
  
  # 1. Check package extdata
  extdata_dir <- system.file("extdata", package = "barcurateR")
  if (nzchar(extdata_dir) && dir.exists(extdata_dir)) {
    ext_files <- list.files(extdata_dir, pattern = "\\.(sqlite|db)$", full.names = TRUE, ignore.case = TRUE)
    paths <- c(paths, ext_files)
  }
  
  # 2. Check barcurateR cache directory
  cache_dir <- tryCatch(rb_db_dir(create = FALSE), error = function(e) "")
  if (nzchar(cache_dir) && dir.exists(cache_dir)) {
    cache_files <- list.files(cache_dir, pattern = "\\.(sqlite|db)$", full.names = TRUE, ignore.case = TRUE)
    paths <- c(paths, cache_files)
  }
  
  # 3. Check legacy regionbarcoder cache directory
  legacy_dir <- tryCatch(tools::R_user_dir("regionbarcoder", which = "data"), error = function(e) "")
  if (nzchar(legacy_dir) && dir.exists(legacy_dir)) {
    legacy_files <- list.files(legacy_dir, pattern = "\\.(sqlite|db)$", full.names = TRUE, ignore.case = TRUE)
    paths <- c(paths, legacy_files)
  }
  
  paths <- unique(paths)
  
  # Filter to only include valid SQLite files
  paths <- paths[vapply(paths, rb_is_sqlite_file, logical(1))]
  
  if (length(paths) == 0) {
    return(data.frame(
      name = character(0),
      version = character(0),
      path = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  names_vec <- tools::file_path_sans_ext(basename(paths))
  versions <- ifelse(grepl("YZFishDB", basename(paths), ignore.case = TRUE), 
                     "full-local", 
                     ifelse(grepl("small_refdb|demo", basename(paths), ignore.case = TRUE),
                            "demo", 
                            "user-curated"))
  
  data.frame(
    name = names_vec,
    version = versions,
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}
