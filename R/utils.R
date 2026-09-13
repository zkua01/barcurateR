rb_match_arg <- function(x, choices, arg = deparse(substitute(x))) {
  if (length(x) != 1 || is.na(x)) {
    stop(arg, " must be a single non-missing value.", call. = FALSE)
  }
  match.arg(tolower(x), choices)
}

rb_quote_in <- function(con, values) {
  if (is.null(values)) return(NULL)
  values <- unique(as.character(values))
  paste(DBI::dbQuoteString(con, values), collapse = ",")
}

rb_clean_sequence <- function(x) {
  toupper(gsub("\\s+", "", as.character(x)))
}

rb_required_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

rb_find_full_yzfishdb <- function() {
  env_path <- Sys.getenv("RB_YZFISHDB_PATH", "")
  if (nzchar(env_path) && file.exists(env_path)) {
    return(normalizePath(env_path, winslash = "/", mustWork = TRUE))
  }
  cached_path <- rb_db_path()
  if (file.exists(cached_path)) {
    return(normalizePath(cached_path, winslash = "/", mustWork = TRUE))
  }
  roots <- normalizePath(c(getwd(), dirname(getwd()), dirname(dirname(getwd())),
                           dirname(dirname(dirname(getwd())))),
                         winslash = "/", mustWork = FALSE)
  candidates <- file.path(roots, "zenodo", "data", "YZFishDB.db")
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) return("")
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

rb_demo_db_path <- function() {
  demo_csv <- system.file("extdata", "yzfishdb_demo.csv", package = "barcurateR")
  if (!nzchar(demo_csv)) {
    demo_csv <- file.path(getwd(), "inst", "extdata", "yzfishdb_demo.csv")
  }
  if (!file.exists(demo_csv)) return("")
  demo_db <- file.path(tempdir(), "barcurateR_demo.sqlite")
  if (!file.exists(demo_db)) {
    demo <- utils::read.csv(demo_csv, stringsAsFactors = FALSE, check.names = FALSE)
    con <- DBI::dbConnect(RSQLite::SQLite(), demo_db)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbWriteTable(con, "yzfishdb_final", demo, overwrite = TRUE)
    DBI::dbWriteTable(con, "yzfishdb_raw", demo[, c("source", "species", "seq_type", "sequence_id", "unique_code", "sequence"), drop = FALSE], overwrite = TRUE)
    qc <- data.frame(
      source = demo$source,
      sequence_id = demo$sequence_id,
      unique_code = demo$unique_code,
      species = demo$species,
      gene = demo$seq_type,
      sequence = demo$sequence,
      qc_flag = demo$qc_flag,
      stringsAsFactors = FALSE
    )
    DBI::dbWriteTable(con, "qc_reference_p2", qc, overwrite = TRUE)
  }
  demo_db
}

rb_find_yzfishdb_data_dir <- function() {
  env_path <- Sys.getenv("RB_YZFISHDB_DATA_DIR", "")
  if (nzchar(env_path) && dir.exists(env_path)) {
    return(normalizePath(env_path, winslash = "/", mustWork = TRUE))
  }
  roots <- normalizePath(c(getwd(), dirname(getwd()), dirname(dirname(getwd())),
                           dirname(dirname(dirname(getwd())))),
                         winslash = "/", mustWork = FALSE)
  candidates <- file.path(roots, "zenodo", "data")
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0) return("")
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

rb_default_data_file <- function(filename) {
  data_dir <- rb_find_yzfishdb_data_dir()
  if (nzchar(data_dir)) {
    candidate <- file.path(data_dir, filename)
    if (file.exists(candidate)) return(candidate)
  }
  bundled <- system.file("extdata", filename, package = "barcurateR")
  if (nzchar(bundled)) return(bundled)
  local <- file.path(getwd(), "inst", "extdata", filename)
  if (file.exists(local)) return(local)
  ""
}

# ------------------------------------------------------------
# rb_standardize_columns()
#
# Renames raw source columns to standardized column names.
#
# column_map should be of the form:
#   c(standard_name = "raw_name")
#
# Example:
#   c(
#     sequence_id = "accession",
#     species = "scientific_name",
#     sequence = "sequence"
#   )
# ------------------------------------------------------------

rb_standardize_columns <- function(raw, column_map) {
  if (is.null(raw) || ncol(raw) == 0) {
    return(raw)
  }
  
  if (is.null(column_map) || length(column_map) == 0) {
    return(raw)
  }
  
  # Allow named lists as well as named character vectors.
  if (is.list(column_map)) {
    column_map <- unlist(column_map, use.names = TRUE)
  }
  
  # Capture names before any further coercion.
  map_names <- names(column_map)
  map_values <- unname(column_map)
  
  if (
    is.null(map_names) ||
    length(map_names) != length(map_values) ||
    any(!nzchar(map_names))
  ) {
    stop(
      "column_map must be a named vector: c(standard_name = \"raw_name\").",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(map_names)) {
    stop("column_map names must be unique.", call. = FALSE)
  }
  
  map_values <- as.character(map_values)
  names(map_values) <- map_names
  
  missing_raw_columns <- setdiff(map_values, names(raw))
  
  if (length(missing_raw_columns) > 0) {
    stop(
      "column_map refers to column(s) not present in the raw data: ",
      paste(missing_raw_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  out <- raw
  
  for (standard_name in map_names) {
    raw_name <- map_values[[standard_name]]
    
    if (!identical(raw_name, standard_name)) {
      out[[standard_name]] <- out[[raw_name]]
    }
  }
  
  # Drop original raw columns that were renamed, unless they already
  # have the standardized name.
  columns_to_drop <- setdiff(map_values, map_names)
  
  out <- out[, setdiff(names(out), columns_to_drop), drop = FALSE]
  
  out
}