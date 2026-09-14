# ============================================================
# data-raw/build-small-demo-data.R
#
# Builds the small demo SQLite database from the committed
# fixture files using the full curation pipeline.
#
#
# Output:
#   inst/extdata/small_refdb.sqlite
# ============================================================

devtools::load_all()

# ----------------------------------------------------------
# 1. Locate fixture files
# ----------------------------------------------------------

ncbi_path <- file.path("tests", "testthat", "fixtures", "ncbi_marker_sequences.csv")
bold_path <- file.path("tests", "testthat", "fixtures", "bold_coi_sequences.tsv")
tax_path  <- file.path("tests", "testthat", "fixtures", "quick_taxonomy.csv")

stopifnot(file.exists(ncbi_path))
stopifnot(file.exists(bold_path))
stopifnot(file.exists(tax_path))

# ----------------------------------------------------------
# 2. Load fixture data
# ----------------------------------------------------------

ncbi_raw <- readr::read_csv(ncbi_path, show_col_types = FALSE)
bold_raw <- readr::read_tsv(bold_path, show_col_types = FALSE)
taxonomy_table <- readr::read_csv(tax_path, show_col_types = FALSE)

message("Loaded ", nrow(ncbi_raw), " NCBI records")
message("Loaded ", nrow(bold_raw), " BOLD records")
message("Loaded ", nrow(taxonomy_table), " taxonomy entries")

# ----------------------------------------------------------
# 3. Clean NCBI fixture
# ----------------------------------------------------------

# Remove predicted / non-reference records
bad_pattern <- paste0(
  "PREDICTED|mRNA|synthetic construct|",
  "uncultured|environmental sample|vector"
)

ncbi_raw <- ncbi_raw[
  !grepl(bad_pattern, ncbi_raw$description, ignore.case = TRUE),
  ,
  drop = FALSE
]

# Keep only one record per accession
ncbi_raw <- ncbi_raw[!duplicated(ncbi_raw$sequence_id), , drop = FALSE]

# Use marker_query as the description for marker extraction
ncbi_raw$description <- ncbi_raw$marker_query

message("After cleaning: ", nrow(ncbi_raw), " NCBI records")

# ----------------------------------------------------------
# 4. Prepare BOLD fixture
# ----------------------------------------------------------

bold_raw$description <- "COI"

# ----------------------------------------------------------
# 5. Parse and standardize sources
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

message("Parsed ", nrow(ncbi_parsed), " NCBI sequences")
message("Parsed ", nrow(bold_parsed), " BOLD sequences")

# ----------------------------------------------------------
# 6. Combine sources
# ----------------------------------------------------------

combined <- rb_combine_sources(
  list(ncbi_parsed, bold_parsed),
  species_col = "species",
  sequence_col = "sequence"
)

message("Combined: ", nrow(combined), " sequences")

# ----------------------------------------------------------
# 7. Resolve ambiguities
# ----------------------------------------------------------

combined <- rb_resolve_ambiguous(
  combined,
  on_unresolved = "drop"
)

message("After ambiguity resolution: ", nrow(combined), " sequences")

# ----------------------------------------------------------
# 8. Run curation pipeline
# ----------------------------------------------------------

out_db <- file.path("inst", "extdata", "small_refdb.sqlite")
dir.create(dirname(out_db), recursive = TRUE, showWarnings = FALSE)

if (file.exists(out_db)) {
  unlink(out_db)
}

result <- rb_curate_reference(
  data = combined,
  taxonomy_table = taxonomy_table,
  run_classifier = FALSE,
  run_divergence = FALSE,
  run_barcode_gap = FALSE,
  db_path = out_db,
  check_ambiguity = TRUE,
  on_ambiguous = "stop",
  require_taxonomy = TRUE
)

# ----------------------------------------------------------
# 9. Verify output
# ----------------------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), out_db)
tables <- DBI::dbListTables(con)
DBI::dbDisconnect(con)

message("\n==================================================")
message("Demo database created: ", out_db)
message("Tables: ", paste(tables, collapse = ", "))
message("Final sequences: ", nrow(result$final_data))
message("==================================================")

# Print species and marker summary
species_summary <- table(result$final_data$species)
marker_summary <- table(result$final_data$seq_type)

message("\nSpecies distribution:")
print(species_summary)

message("\nMarker distribution:")
print(marker_summary)