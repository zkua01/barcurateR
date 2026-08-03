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
  demo_csv <- system.file("extdata", "yzfishdb_demo.csv", package = "regionbarcoder")
  if (!nzchar(demo_csv)) {
    demo_csv <- file.path(getwd(), "inst", "extdata", "yzfishdb_demo.csv")
  }
  if (!file.exists(demo_csv)) return("")
  demo_db <- file.path(tempdir(), "regionbarcoder_demo.sqlite")
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
  bundled <- system.file("extdata", filename, package = "regionbarcoder")
  if (nzchar(bundled)) return(bundled)
  local <- file.path(getwd(), "inst", "extdata", filename)
  if (file.exists(local)) return(local)
  ""
}





# ============================================================
# MODULE 5 — Taxonomy join, QC filtering, column standardization
# Destination guide is noted above each function/fix.
# ============================================================
 
#' Standardize column names to the package's internal schema
#'
#' Renames arbitrary source-table column names to the standard names the
#' rest of regionbarcoder expects (species, sequence, qc_flag, etc.). This
#' is the piece that lets Module 1's per-source parsing (NCBI/BOLD/MitoFish/
#' MIDORI2/local, or any future source) hand off to the rest of the pipeline
#' without every downstream function needing to know each source's original
#' column names.
#'
#' @param data A data.frame with arbitrary column names.
#' @param column_map Named character vector: names(column_map) = the
#'   standard name expected internally, values = the column name actually
#'   present in `data`. e.g.
#'   c(species = "species_name", sequence = "seq", qc_flag = "flag")
#' @param drop_unmapped If FALSE (default), columns not mentioned in
#'   column_map are kept as-is. If TRUE, only the mapped columns are kept
#'   (useful right before writing to a shared table like yzfishdb_final).
#' @return `data` with columns renamed to the standard schema.
rb_standardize_columns <- function(data, column_map, drop_unmapped = FALSE) {
  missing_source <- setdiff(unname(column_map), names(data))
  if (length(missing_source) > 0) {
    stop("Columns referenced in column_map not found in data: ",
         paste(missing_source, collapse = ", "), call. = FALSE)
  }
 
  # column_map is standard_name = source_name; we rename source_name -> standard_name
  for (standard_name in names(column_map)) {
    source_name <- column_map[[standard_name]]
    if (!identical(source_name, standard_name)) {
      names(data)[names(data) == source_name] <- standard_name
    }
  }
 
  if (isTRUE(drop_unmapped)) {
    data <- data[, names(column_map), drop = FALSE]
  }
  data
}
 
# --- Test ---
# raw <- data.frame(species_name = "Homo sapiens", seq = "ACGT", flag = "pass")
# rb_standardize_columns(raw, column_map = c(species = "species_name",
#                                            sequence = "seq",
#                                            qc_flag = "flag"))
# expect: columns renamed to species / sequence / qc_flag; values untouched
 
