#' Install and manage barcurateR reference databases
#'
#' @description
#' Functions to locate, download, and manage cached reference databases.
#' 
#' * `rb_db_dir()`: Returns the path to the local cache directory.
#' * `rb_db_path()`: Returns the full path to a specific database file in the cache.
#' * `rb_db_available()`: Checks if a database file exists locally.
#' * `rb_db_url()`: Returns the configured download URL for the full database.
#' * `rb_install_db()`: Downloads the full database from Zenodo (or a custom URL) 
#'   into the local cache, with retry logic and checksum verification.
#' * `rb_yzfishdb_release()`: Returns metadata for the official YZFishDB Zenodo release.
#' * `rb_default_yzfishdb_url()`: Returns the default Zenodo API URL for YZFishDB.
#'
#' @param create Logical. If `TRUE`, create the cache directory if it doesn't exist.
#' @param filename The name of the database file (default `"YZFishDB.db"`).
#' @param path Path to check for database availability.
#' @param url URL to download the database from.
#' @param destfile Destination file path for the download.
#' @param overwrite Logical. Overwrite existing file?
#' @param md5,sha256 Expected checksums for verification.
#' @param min_bytes Minimum expected file size in bytes.
#' @param quiet Logical. Suppress download progress messages?
#' @param timeout Download timeout in seconds (default 3600).
#' @param retries Number of download retry attempts.
#'
#' @return 
#' * `rb_db_dir()` and `rb_db_path()` return a character string path.
#' * `rb_db_available()` returns a logical value.
#' * `rb_install_db()` returns the normalized path to the downloaded file.
#' * `rb_yzfishdb_release()` returns a data frame of release metadata.
#'
#' @details
#' The cache directory is determined by:
#' 1. The `BARCURATER_CACHE_DIR` environment variable.
#' 2. The default R user data directory (`tools::R_user_dir("barcurateR", "data")`).
#'
#' The download URL is determined by:
#' 1. The `BARCURATER_DB_URL` environment variable.
#' 2. The `barcurateR.db_url` global option.
#' 3. The default YZFishDB Zenodo URL.
#'
#' @name rb_db_install
#' @family database management
NULL

#' @rdname rb_db_install
#' @export
rb_yzfishdb_release <- function() {
  data.frame(
    doi = "10.5281/zenodo.18155084",
    record = "18155084",
    filename = "YZFishDB.db",
    url = rb_default_yzfishdb_url(),
    size_bytes = 635740160,
    md5 = "0e5e0c3e294c4a55bcc3065c26db9c84",
    stringsAsFactors = FALSE
  )
}

#' @rdname rb_db_install
#' @export
rb_default_yzfishdb_url <- function() {
  "https://zenodo.org/api/records/18155084/files/YZFishDB.db/content"
}

#' @rdname rb_db_install
#' @export
rb_db_url <- function() {
  env_url <- Sys.getenv("BARCURATER_DB_URL", "")
  if (nzchar(env_url)) return(env_url)
  
  opt_url <- getOption("barcurateR.db_url", rb_default_yzfishdb_url())
  if (is.null(opt_url)) return("")
  
  as.character(opt_url)
}

#' @rdname rb_db_install
#' @export
rb_db_dir <- function(create = TRUE) {
  cache_dir <- Sys.getenv("BARCURATER_CACHE_DIR", "")
  
  if (!nzchar(cache_dir)) {
    cache_dir <- tools::R_user_dir("barcurateR", which = "data")
  }
  
  if (isTRUE(create) && !dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  normalizePath(cache_dir, winslash = "/", mustWork = dir.exists(cache_dir))
}

#' @rdname rb_db_install
#' @export
rb_db_path <- function(filename = "YZFishDB.db") {
  file.path(rb_db_dir(), filename)
}

#' @rdname rb_db_install
#' @export
rb_db_available <- function(path = rb_db_path()) {
  nzchar(path) && file.exists(path)
}

#' @rdname rb_db_install
#' @export
rb_install_db <- function(url = rb_db_url(), 
                          destfile = rb_db_path(),
                          overwrite = FALSE, 
                          md5 = NULL, 
                          sha256 = NULL,
                          min_bytes = 1, 
                          quiet = FALSE,
                          timeout = 3600, 
                          retries = 3) {
  
  if (file.exists(destfile) && !isTRUE(overwrite)) {
    return(normalizePath(destfile, winslash = "/", mustWork = TRUE))
  }
  
  if (!nzchar(url)) {
    stop(
      "No database download URL is configured. ",
      "Set BARCURATER_DB_URL, set option barcurateR.db_url, ",
      "or pass url = '...' to rb_install_db().",
      call. = FALSE
    )
  }
  
  dest_dir <- dirname(destfile)
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  tmp <- tempfile("barcurateR_db_", fileext = ".db")
  on.exit(unlink(tmp), add = TRUE)
  
  rb_download_with_retries(
    url = url,
    destfile = tmp,
    quiet = quiet,
    timeout = timeout,
    retries = retries,
    downloader = getOption("barcurateR.download_file", rb_download_file)
  )
  
  # Auto-fill checksums/size if downloading the official YZFishDB release
  if (identical(url, rb_default_yzfishdb_url())) {
    release <- rb_yzfishdb_release()
    if (is.null(md5)) md5 <- release$md5[[1]]
    if (identical(min_bytes, 1)) min_bytes <- release$size_bytes[[1]]
  }
  
  file_size <- file.info(tmp)$size
  if (is.na(file_size) || file_size < min_bytes) {
    stop("Downloaded database file is empty or incomplete.", call. = FALSE)
  }
  
  if (!is.null(md5)) {
    observed <- unname(tools::md5sum(tmp))
    if (!identical(tolower(observed), tolower(md5))) {
      stop("Downloaded database md5 checksum does not match expected md5.", call. = FALSE)
    }
  }
  
  if (!is.null(sha256)) {
    observed <- unname(tools::sha256sum(tmp))
    if (!identical(tolower(observed), tolower(sha256))) {
      stop("Downloaded database checksum does not match expected sha256.", call. = FALSE)
    }
  }
  
  ok <- file.copy(tmp, destfile, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Could not write database to: ", destfile, call. = FALSE)
  }
  
  normalizePath(destfile, winslash = "/", mustWork = TRUE)
}

#' Internal download helper with retries
#' @noRd
rb_download_with_retries <- function(url, destfile, quiet, timeout, retries, downloader) {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(as.numeric(timeout), as.numeric(old_timeout), na.rm = TRUE))
  
  retries <- max(1L, as.integer(retries))
  last_error <- NULL
  
  for (attempt in seq_len(retries)) {
    unlink(destfile)
    
    result <- tryCatch(
      {
        downloader(url, destfile, quiet = quiet)
        TRUE
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        FALSE
      }
    )
    
    if (isTRUE(result) && file.exists(destfile)) {
      return(invisible(destfile))
    }
    
    if (!quiet && attempt < retries) {
      message("Download failed; retrying (", attempt + 1L, "/", retries, ").")
    }
  }
  
  stop(
    "Could not download database after ", retries, " attempt(s). ",
    "For large files, try rb_install_db(timeout = 7200) or download manually. ",
    "Last error: ", last_error,
    call. = FALSE
  )
}

#' Internal file/URL download handler
#' @noRd
rb_download_file <- function(url, destfile, quiet = FALSE) {
  # Handle local file paths
  if (file.exists(url)) {
    ok <- file.copy(url, destfile, overwrite = TRUE)
    if (!isTRUE(ok)) stop("Could not copy local file: ", url, call. = FALSE)
    return(invisible(destfile))
  }
  
  # Handle file:// URLs
  if (grepl("^file://", url)) {
    local_path <- utils::URLdecode(sub("^file:///?", "", url))
    if (.Platform$OS.type == "windows") {
      local_path <- sub("^/([A-Za-z]):", "\\1:", local_path)
    }
    ok <- file.copy(local_path, destfile, overwrite = TRUE)
    if (!isTRUE(ok)) stop("Could not copy local file: ", local_path, call. = FALSE)
    return(invisible(destfile))
  }
  
  # Handle HTTP/HTTPS URLs
  status <- utils::download.file(url, destfile, mode = "wb", quiet = quiet)
  if (!identical(status, 0L)) {
    stop("Could not download from URL: ", url, call. = FALSE)
  }
  
  invisible(destfile)
}