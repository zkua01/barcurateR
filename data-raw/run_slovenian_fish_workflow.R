# ============================================================
# Slovenian Freshwater Fish eDNA Reference Database Workflow
# ============================================================
# This script demonstrates the full barcurateR pipeline:
# 1. Standardization of heterogeneous NCBI and BOLD downloads
# 2. Ambiguity detection and resolution
# 3. The 5-stage curation pipeline (QC, ML, Taxonomy, Gap, Export)
# 4. Database auditing and diagnostics
# 5. Export for downstream metabarcoding pipelines
# 6. Ecological validation via eDNA taxonomic assignment
# ============================================================

library(barcurateR)
library(dplyr)
library(readr)
library(tidyr)

# Create directories for outputs
dir.create("exports", showWarnings = FALSE)
dir.create("data-raw/slovenian_fish/processed", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# STAGE 0: Load Raw Data and Taxonomy
# ----------------------------------------------------------
message(">>> Loading raw data...")

# Load the raw downloads (assuming they were saved as CSVs by your download script)
ncbi_raw <- read_csv("data-raw/slovenian_fish/raw/ncbi_sequence_raw.csv", show_col_types = FALSE) |>
  filter(!is.na(sequence))
bold_raw <- read_csv("data-raw/slovenian_fish/raw/bold_sequence_raw.csv", show_col_types = FALSE) |>
  filter(!is.na(sequence))

# Load the curated taxonomy table for the 99 Slovenian fish species
taxonomy_table <- read_csv("data-raw/slovenian_fish/processed/slovenian_fish_taxonomy.csv", show_col_types = FALSE)

message(sprintf("Loaded %d NCBI records and %d BOLD records.", nrow(ncbi_raw), nrow(bold_raw)))


# ----------------------------------------------------------
# STAGE 1: Source Standardization
# ----------------------------------------------------------
message("\n>>> Standardizing source tables...")

# 1A. NCBI Standardization
# We extract the accession from the header to use as a clean sequence_id,
# and pass the full header to description_col so rb_extract_markers() 
# can parse out "12S", "16S", "COI", or "complete_genome".
ncbi_raw$sequence_id <- sub("\\s.*", "", ncbi_raw$header)

ncbi_parsed <- rb_parse_source_table(
  raw = ncbi_raw,
  source_name = "ncbi",
  column_map = c(
    sequence_id = "sequence_id",
    species = "species",
    sequence = "sequence"
  ),
  description_col = "header"
)

# 1B. BOLD Standardization
# BOLD uses "COI-5P". We pass the seq_type column to description_col 
# so rb_simplify_seq_type() converts it to standard "COI".
bold_raw$description <- bold_raw$seq_type 

bold_parsed <- rb_parse_source_table(
  raw = bold_raw,
  source_name = "bold",
  column_map = c(
    sequence_id = "sequence_id",
    species = "species_name",
    sequence = "sequence"
  ),
  description_col = "description"
)

message(sprintf("Standardized: %d NCBI seqs, %d BOLD seqs.", nrow(ncbi_parsed), nrow(bold_parsed)))


# ----------------------------------------------------------
# STAGE 2: Integration & Ambiguity Resolution
# ----------------------------------------------------------
message("\n>>> Combining sources and resolving ambiguities...")

combined <- rb_combine_sources(list(ncbi_parsed, bold_parsed))

# Detect sequences assigned to multiple species
combined <- rb_detect_ambiguity(combined)
n_ambig <- sum(combined$is_ambiguous)
message(sprintf("Found %d ambiguous sequences.", n_ambig))

if (n_ambig > 0) {
  # Archive ambiguities for manual review/auditing
  rb_write_ambiguity_csv(
    combined[combined$is_ambiguous, ], 
    "data-raw/slovenian_fish/processed/ambiguous_sequences.csv"
  )
}

# For this automated run, we drop unresolved ambiguities to prevent 
# downstream QC and ML from failing on conflicting labels.
resolved <- rb_resolve_ambiguous(combined, on_unresolved = "drop")


# ----------------------------------------------------------
# STAGE 3: The Curation Pipeline (QC, ML, Taxonomy, SQLite)
# ----------------------------------------------------------
message("\n>>> Running the main curation pipeline...")

db_path <- file.path("data-raw/slovenian_fish/processed/slovenian_fish_refdb.sqlite")
contam_path <- file.path("tests/testthat/fixtures/contaminants/contaminant_db")
numts_path <- file.path("tests/testthat/fixtures/contaminants/fish_numts.fasta")

# Note: blast_db and numt_fasta are set to NULL here for the demo.
# In a production run, you would point these to your contaminant/NUMT databases.
curated <- rb_curate_reference(
  data = resolved,
  blast_db = contam_path,
  numt_fasta = numts_path,
  run_classifier = TRUE,       # Train ML to reclassify "other" sequences
  run_divergence = TRUE,      # Requires mafft/FastTree (skip for speed)
  run_barcode_gap = TRUE,      # Calculate genetic distances & diagnostic sites
  taxonomy_table = taxonomy_table,
  db_path = db_path,
  on_ambiguous = "stop",
  require_taxonomy = TRUE,
  archive_ambiguity = TRUE,
  ambiguity_path = "data-raw/slovenian_fish/processed/ambiguous_sequences.csv"
)

message(sprintf("Curation complete. Final database contains %d sequences.", nrow(curated$final_data)))


# ----------------------------------------------------------
# STAGE 4: Database Auditing & Diagnostics
# ----------------------------------------------------------
message("\n>>> Auditing the curated database...")

con <- rb_connect(db_path)

# 4A. Overall QC Health
qc_sum <- rb_qc_summary(con)
print("QC Flag Distribution:")
print(qc_sum)

# 4B. Marker Coverage by Taxonomic Rank
marker_cov <- rb_marker_coverage(con, rank = "order")
print("Marker Coverage by Order:")
print(head(marker_cov))

# 4C. Taxonomic Coverage
taxa_cov <- rb_list_taxa(con, rank = "species")
print("Taxonomic Coverage by Species:")
print(head(taxa_cov))

# 4C. Barcode Gap Metrics (Identify species with poor genetic resolution)
gap_metrics <- rb_barcode_gap(con)
if (nrow(gap_metrics) > 0) {
  poor_resolution <- gap_metrics %>% 
    filter(taxonomic_resolution != "species" | assignment_risk == "high")
  message(sprintf("Warning: %d species/markers have high assignment risk or poor resolution.", nrow(poor_resolution)))
}

head(poor_resolution)

rb_disconnect(con)

# export audit data frames for vignette
audit_output <- list(qc_sum, marker_cov, taxa_cov, poor_resolution)
saveRDS(audit_output, "inst/extdata/case_study_outputs.rds")

# ----------------------------------------------------------
# STAGE 5: Downstream Export
# ----------------------------------------------------------
message("\n>>> Exporting for downstream metabarcoding pipelines...")

con <- rb_connect(db_path)

# Extract only high-quality, passing sequences
refs_12S <- rb_get_sequences(con, marker = "12S", qc_flag = "pass")
refs_COI <- rb_get_sequences(con, marker = "COI", qc_flag = "pass")

# Export 12S for DADA2 (e.g., MiFish/Teleo primers)
rb_export_dada2(refs_12S, "exports/slovenian_fish_12S_dada2.fasta")

# Export COI for QIIME2 (generates both FASTA and Taxonomy TSV)
rb_export_qiime2(
  refs_COI, 
  "exports/slovenian_fish_COI_qiime2.fasta", 
  "exports/slovenian_fish_COI_taxonomy.tsv"
)

# Export 12S database for BLASTn
rb_export_blastn(
  rb_get_sequences(con, marker = "12S", qc_flag = "pass"), 
  "exports/slovenian_fish_12s_blastn.fasta"
)

# Export full database for BLASTn
rb_export_blastn(
  rb_get_sequences(con, qc_flag = "pass"), 
  "exports/slovenian_fish_all_markers_blastn.fasta"
)

rb_disconnect(con)
message("Exports saved to the 'exports/' directory.")


# ----------------------------------------------------------
# STAGE 6: Ecological Validation (Sava River eDNA Assignment)
# ----------------------------------------------------------

message("\n>>> Validating against Sava River eDNA data...")

sava_asvs_fasta <- "data-raw/sava_river/processed/sava_asvs.fasta"
sava_counts <- read_csv("data-raw/sava_river/processed/sava_abundance_table.csv", show_col_types = FALSE)

con <- rb_connect(db_path)

# Run BLASTn assignment against the curated 12S database (1st pass)
res_pass1 <- rb_assign_edna(
  asv_fasta = sava_asvs_fasta,
  con = con,
  marker = "12S",
  method = "blastn",
  min_identity = 99,
  out_dir = "exports/blastn_pass1"
)

all_asv_ids    <- names(Biostrings::readDNAStringSet(sava_asvs_fasta))
assigned_ids   <- res_pass1$asv_id[res_pass1$assignment_status != "no_match"]
unassigned_ids <- setdiff(all_asv_ids, assigned_ids)

message(sprintf("Pass 1: %d assigned, %d unassigned -> pass 2.",
                length(assigned_ids), length(unassigned_ids)))
if (length(unassigned_ids) > 0) {
  all_seqs <- Biostrings::readDNAStringSet(sava_asvs_fasta)
  pass2_fasta <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(all_seqs[names(all_seqs) %in% unassigned_ids], pass2_fasta)
  res_pass2 <- rb_assign_edna(
    asv_fasta = pass2_fasta, con = con, marker = NULL,
    method = "blastn", qc_flag = c("pass","manual_review_needed"), min_identity = 99, min_coverage = 0.9,
    max_target_seqs = 500, out_dir = "exports/blastn_pass2"
  )
  res_pass2$db_source <- "full_reference"
} else {
  res_pass2 <- res_pass1[0, ]
}

final_assignments <- bind_rows(
  res_pass1 %>% filter(assignment_status != "no_match") %>%
    mutate(db_source = "12S_subset"),
  res_pass2
)

missing_ids <- setdiff(all_asv_ids, final_assignments$asv_id)
if (length(missing_ids) > 0) {
  pad <- final_assignments[rep(NA_integer_, length(missing_ids)), , drop = FALSE]
  pad$asv_id            <- missing_ids
  pad$assignment_status <- "no_match"
  pad$db_source         <- "no_hit"
  final_assignments <- bind_rows(final_assignments, pad)
}
stopifnot(nrow(final_assignments) == length(all_asv_ids),
          !any(duplicated(final_assignments$asv_id)))

rb_disconnect(con)

message("Final assignment summary:")
print(table(final_assignments$assignment_status))
print(table(final_assignments$db_source))

write_csv(final_assignments, "inst/extdata/sava_river_12s_assignments_two_pass.csv")

# Build the final Species-by-Sample ecological matrix
community_matrix <- rb_build_species_matrix(
  assignment = final_assignments,
  asv_table = sava_counts,
  include_ambiguous = FALSE,
  unassigned = FALSE
)

# Save the ecological matrix for downstream ordination (e.g., NMDS, PCoA)
write_csv(community_matrix, "inst/extdata/sava_river_community_matrix.csv")

rb_disconnect(con)
message("Ecological validation complete.")
