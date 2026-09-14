#' Detect and archive ambiguous sequences
#'
#' @description
#' Utility functions for detecting sequences assigned to multiple species
#' and archiving ambiguity records for downstream auditing.
#'
#' * `rb_detect_ambiguity()`: Flags sequences that are assigned to more
#'   than one species.
#' * `rb_write_ambiguity_csv()`: Writes an ambiguity table to CSV.
#' * `rb_read_ambiguity_csv()`: Reads an ambiguity table from CSV.
#'
#' @param data A standardized data frame containing at least the columns
#'   specified by `species_col` and `sequence_col`.
#' @param species_col Name of the species column. Default `"species"`.
#' @param sequence_col Name of the sequence column. Default `"sequence"`.
#' @param path Path to the ambiguity CSV file. Default
#'   `"ambiguous_sequences.csv"`.
#'
#' @return
#' * `rb_detect_ambiguity()` returns the input data frame with an
#'   `is_ambiguous` logical column added.
#' * `rb_write_ambiguity_csv()` returns the output path invisibly.
#' * `rb_read_ambiguity_csv()` returns a data frame.
#'
#' @details
#' Ambiguous sequences are sequences that occur under more than one
#' species name. These should normally be resolved before QC using
#' `rb_resolve_ambiguous()` or a manual resolution table.
#'
#' The ambiguity table is an auditing artifact. It is not used directly
#' by the QC pipeline, but it can be archived in the final SQLite
#' database by `rb_curate_reference()` using the `ambiguity_data` or
#' `ambiguity_path` arguments.
#'
#' @examples
#' \dontrun{
#' combined <- rb_combine_sources(list(ncbi_parsed, bold_parsed))
#'
#' combined <- rb_detect_ambiguity(combined)
#'
#' ambiguities <- combined[combined$is_ambiguous, ]
#'
#' if (nrow(ambiguities) > 0) {
#'   rb_write_ambiguity_csv(ambiguities, "ambiguous_sequences.csv")
#' }
#' }
#'
#' @name rb_ambiguity_utils
#' @family parsing
NULL

#' Match a single argument to a set of choices
#'
#' Internal helper for validating single-choice character arguments.
#'
#' @keywords internal
#' @noRd
rb_match_arg <- function(x, choices, arg = deparse(substitute(x))) {
  if (length(x) != 1 || is.na(x)) {
    stop(arg, " must be a single non-missing value.", call. = FALSE)
  }
  match.arg(tolower(x), choices)
}

#' Quote values for use in SQL IN clauses
#'
#' Internal helper for building safe SQL IN conditions.
#'
#' @keywords internal
#' @noRd
rb_quote_in <- function(con, values) {
  if (is.null(values)) return(NULL)
  values <- unique(as.character(values))
  paste(DBI::dbQuoteString(con, values), collapse = ",")
}

#' Clean a DNA sequence string
#'
#' Internal helper that removes whitespace and converts sequences to
#' uppercase.
#'
#' @keywords internal
#' @noRd
rb_clean_sequence <- function(x) {
  toupper(gsub("\\s+", "", as.character(x)))
}

#' Check that required columns are present
#'
#' Internal helper that stops with an informative error if required
#' columns are missing from a data frame.
#'
#' @keywords internal
#' @noRd
rb_required_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Locate the full local YZFishDB database
#'
#' Internal helper that searches common locations for the full
#' YZFishDB SQLite database.
#'
#' @keywords internal
#' @noRd
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

#' Build a temporary demo database from bundled demo data
#'
#' Internal helper for creating a small demo SQLite database.
#'
#' @keywords internal
#' @noRd
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

#' Locate the local YZFishDB data directory
#'
#' Internal helper that searches common locations for YZFishDB
#' supplementary data files.
#'
#' @keywords internal
#' @noRd
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

#' Locate a default bundled data file
#'
#' Internal helper that searches the YZFishDB data directory,
#' installed package `extdata`, and local `inst/extdata` for a file.
#'
#' @keywords internal
#' @noRd
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

#' Rename raw source columns to standardized names
#'
#' Internal helper used by `rb_parse_source_table()` to map raw
#' user-supplied column names onto the standardized barcurateR schema.
#'
#' @keywords internal
#' @noRd
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

#' Check that required taxonomy columns are present
#'
#' Internal helper that enforces the kingdom-to-species taxonomy
#' requirement for final reference tables.
#'
#' @keywords internal
#' @noRd
rb_required_taxonomy_columns <- function(data) {
  required <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(
      "Missing required taxonomy columns: ",
      paste(missing, collapse = ", "),
      ". The final reference table requires full taxonomy from species to kingdom. ",
      "Only occurrence and habitat are optional.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

#' @rdname rb_ambiguity_utils
#' @export
rb_detect_ambiguity <- function(data,
                                species_col = "species",
                                sequence_col = "sequence") {
  rb_required_columns(data, c(species_col, sequence_col))
  
  if (nrow(data) == 0) {
    data$is_ambiguous <- logical(0)
    return(data)
  }
  
  seq_clean <- rb_clean_sequence(data[[sequence_col]])
  
  species_by_seq <- split(data[[species_col]], seq_clean)
  
  ambiguous_sequences <- names(species_by_seq)[
    vapply(
      species_by_seq,
      function(x) length(unique(stats::na.omit(x))) > 1,
      logical(1)
    )
  ]
  
  data$is_ambiguous <- seq_clean %in% ambiguous_sequences
  
  data
}

#' @rdname rb_ambiguity_utils
#' @export
rb_write_ambiguity_csv <- function(data,
                                   path = "ambiguous_sequences.csv") {
  utils::write.csv(data, path, row.names = FALSE)
  invisible(path)
}

#' @rdname rb_ambiguity_utils
#' @export
rb_read_ambiguity_csv <- function(path = "ambiguous_sequences.csv") {
  if (!file.exists(path)) {
    stop("Ambiguity file not found: ", path, call. = FALSE)
  }
  
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

