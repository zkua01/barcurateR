# tests/testthat/helper-mock-db.R
#
# Creates an in-memory SQLite database with generic table names
# and realistic columns for testing query and diagnostic functions.

create_mock_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  
  ref_data <- data.frame(
    sequence_id = paste0("ID_", 1:6),
    unique_code = paste0("seq_", 1:6),
    species = c("Species A", "Species A", "Species B", "Species B", "Species C", "Species C"),
    genus = c("GenusA", "GenusA", "GenusB", "GenusB", "GenusC", "GenusC"),
    family = c("Fam1", "Fam1", "Fam1", "Fam1", "Fam2", "Fam2"),
    order = "Ord1",
    class = "Class1",
    phylum = "Phylum1",
    kingdom = "Kingdom1",
    seq_type = c("COI", "12S", "COI", "COI", "12S", "16S"),
    sequence = c("ACGT", "ACGT", "TTTT", "GGGG", "CCCC", "AAAA"),
    source = c("ncbi", "bold", "ncbi", "ncbi", "bold", "bold"),
    qc_flag = c("pass", "pass", "pass", "contaminant", "pass", "pass"),
    occurrence = c("native", "native", "introduced", NA, "native", "native"),
    habitat = c("freshwater", "freshwater", "marine", NA, "freshwater", "freshwater"),
    stringsAsFactors = FALSE
  )
  
  qc_data <- data.frame(
    unique_code = paste0("seq_", 1:6),
    species = c("Species A", "Species A", "Species B", "Species B", "Species C", "Species C"),
    qc_flag = c("pass", "pass", "pass", "contaminant", "pass", "pass"),
    stringsAsFactors = FALSE
  )
  
  gap_data <- data.frame(
    species = c("Species A", "Species B"),
    marker = c("COI", "COI"),
    gap_exists = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  
  amb_data <- data.frame(
    sequence = c("ACGT", "ACGT"),
    species = c("Species A", "Species X"),
    note = c("test ambiguity", "test ambiguity"),
    stringsAsFactors = FALSE
  )
  
  DBI::dbWriteTable(con, "reference_final", ref_data)
  DBI::dbWriteTable(con, "qc_reference", qc_data)
  DBI::dbWriteTable(con, "barcode_gap_metrics", gap_data)
  DBI::dbWriteTable(con, "ambiguous_sequences", amb_data)
  
  con
}