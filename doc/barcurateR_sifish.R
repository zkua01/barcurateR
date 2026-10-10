## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>", eval=FALSE)
library(barcurateR)
library(dplyr)
library(readr)
library(tidyr)
library(Biostrings)
library(RefManageR)

BibOptions(check.entries = FALSE, style = "html", hyperlink = FALSE, cite.style = 'authoryear', max.names = 1, longnamesfirst = F)
myBib <- ReadBib(system.file("extdata", "barcurateR.bib", package = "barcurateR"), check = FALSE)

outputs <- readRDS(system.file("extdata", "case_study_outputs.rds", package = "barcurateR"))
final_assignments <- read_csv(system.file("extdata", "sava_river_12s_assignments_two_pass.csv", package = "barcurateR"))
community_matrix <- read_csv(system.file("extdata", "sava_river_community_matrix.csv", package = "barcurateR"))

## ----message=FALSE, warning=FALSE---------------------------------------------
# library(barcurateR)
# library(dplyr)
# library(readr)
# library(tidyr)
# library(Biostrings)

## -----------------------------------------------------------------------------
# ncbi_csv <- "data-raw/slovenian_fish/raw/ncbi_sequence_raw.csv"
# bold_csv <- "data-raw/slovenian_fish/raw/bold_sequence_raw.csv"
# taxonomy_csv <- "data-raw/slovenian_fish/processed/slovenian_fish_taxonomy.csv"
# 
# ncbi_raw <- read_csv(ncbi_csv, show_col_types = FALSE) |>
#   filter(!is.na(sequence))
# 
# bold_raw <- read_csv(bold_csv, show_col_types = FALSE) |>
#   filter(!is.na(sequence))
# 
# taxonomy_table <- read_csv(taxonomy_csv, show_col_types = FALSE)

## -----------------------------------------------------------------------------
# # NCBI
# ncbi_raw$sequence_id <- sub("\\s.*", "", ncbi_raw$header)
# 
# ncbi_parsed <- rb_parse_source_table(
#   raw = ncbi_raw,
#   source_name = "ncbi",
#   column_map = c(
#     sequence_id = "sequence_id",
#     species = "species",
#     sequence = "sequence"
#   ),
#   description_col = "header"
# )
# 
# # BOLD
# bold_raw$description <- bold_raw$seq_type
# 
# bold_parsed <- rb_parse_source_table(
#   raw = bold_raw,
#   source_name = "bold",
#   column_map = c(
#     sequence_id = "sequence_id",
#     species = "species_name",
#     sequence = "sequence"
#   ),
#   description_col = "description"
# )

## -----------------------------------------------------------------------------
# combined <- rb_combine_sources(list(ncbi_parsed, bold_parsed))

## -----------------------------------------------------------------------------
# combined <- rb_detect_ambiguity(combined)
# n_ambig <- sum(combined$is_ambiguous)
# 
# if (n_ambig > 0) {
#   rb_write_ambiguity_csv(
#     combined[combined$is_ambiguous, ],
#     "data-raw/slovenian_fish/processed/ambiguous_sequences.csv"
#   )
# }
# 
# resolved <- rb_resolve_ambiguous(
#   combined,
#   on_unresolved = "drop"
# )

## -----------------------------------------------------------------------------
# contam_path <- "tests/testthat/fixtures/contaminants/contaminant_db"
# numts_path  <- "tests/testthat/fixtures/contaminants/fish_numts.fasta"
# db_path <- "data-raw/slovenian_fish/processed/slovenian_fish_refdb.sqlite"
# 
# curated <- rb_curate_reference(
#   data = resolved,
#   blast_db = contam_path,
#   numt_fasta = numts_path,
#   run_classifier = TRUE,
#   run_divergence = TRUE,  # requires MAFFT and FastTree
#   run_barcode_gap = TRUE,
#   taxonomy_table = taxonomy_table,
#   db_path = db_path,
#   on_ambiguous = "stop",
#   require_taxonomy = TRUE,
#   archive_ambiguity = TRUE,
#   ambiguity_path = "data-raw/slovenian_fish/processed/ambiguous_sequences.csv"
# )

## -----------------------------------------------------------------------------
# con <- rb_connect(db_path)

## -----------------------------------------------------------------------------
# qc_sum <- rb_qc_summary(con)
# print(qc_sum)

## ----eval=TRUE----------------------------------------------------------------
print(outputs[[1]])

## -----------------------------------------------------------------------------
# marker_cov <- rb_marker_coverage(con, rank = "order")
# head(marker_cov)

## ----eval=TRUE----------------------------------------------------------------
head(outputs[[2]])

## -----------------------------------------------------------------------------
# taxa_cov <- rb_list_taxa(con, rank = "species")
# head(taxa_cov)

## ----eval=TRUE----------------------------------------------------------------
head(outputs[[3]])

## -----------------------------------------------------------------------------
# gap_metrics <- rb_barcode_gap(con)
# 
# if (nrow(gap_metrics) > 0) {
#   poor_resolution <- gap_metrics |>
#     filter(taxonomic_resolution != "species" | assignment_risk == "high")
# }
# 
# head(poor_resolution)

## ----eval=TRUE----------------------------------------------------------------
head(outputs[[4]])

## -----------------------------------------------------------------------------
# rb_disconnect(con)

## -----------------------------------------------------------------------------
# con <- rb_connect(db_path)
# 
# refs_12S <- rb_get_sequences(con, marker = "12S", qc_flag = "pass")
# refs_COI <- rb_get_sequences(con, marker = "COI", qc_flag = "pass")
# 
# # DADA2
# rb_export_dada2(
#   refs_12S,
#   "exports/slovenian_fish_12S_dada2.fasta"
# )
# 
# # QIIME2
# rb_export_qiime2(
#   refs_COI,
#   "exports/slovenian_fish_COI_qiime2.fasta",
#   "exports/slovenian_fish_COI_taxonomy.tsv"
# )
# 
# # BLASTn: 12S only
# rb_export_blastn(
#   rb_get_sequences(con, marker = "12S", qc_flag = "pass"),
#   "exports/slovenian_fish_12s_blastn.fasta"
# )
# 
# # BLASTn: all markers
# rb_export_blastn(
#   rb_get_sequences(con, qc_flag = "pass"),
#   "exports/slovenian_fish_all_markers_blastn.fasta"
# )
# 
# rb_disconnect(con)

## -----------------------------------------------------------------------------
# sava_asvs_fasta <- "data-raw/sava_river/processed/sava_asvs.fasta"
# sava_counts <- read_csv(
#   "data-raw/sava_river/processed/sava_abundance_table.csv",
#   show_col_types = FALSE
# )

## -----------------------------------------------------------------------------
# con <- rb_connect(db_path)
# 
# # Pass 1: 12S subset
# res_pass1 <- rb_assign_edna(
#   asv_fasta = sava_asvs_fasta,
#   con = con,
#   marker = "12S",
#   method = "blastn",
#   min_identity = 99,
#   out_dir = "exports/blastn_pass1"
# )
# 
# all_asv_ids <- names(Biostrings::readDNAStringSet(sava_asvs_fasta))
# assigned_ids <- res_pass1$asv_id[res_pass1$assignment_status != "no_match"]
# unassigned_ids <- setdiff(all_asv_ids, assigned_ids)
# 
# # Pass 2: full reference for unassigned ASVs
# if (length(unassigned_ids) > 0) {
#   all_seqs <- Biostrings::readDNAStringSet(sava_asvs_fasta)
#   pass2_fasta <- tempfile(fileext = ".fasta")
# 
#   Biostrings::writeXStringSet(
#     all_seqs[names(all_seqs) %in% unassigned_ids],
#     pass2_fasta
#   )
# 
#   res_pass2 <- rb_assign_edna(
#     asv_fasta = pass2_fasta,
#     con = con,
#     marker = NULL,
#     method = "blastn",
#     qc_flag = c("pass", "manual_review_needed"),
#     min_identity = 99,
#     min_coverage = 0.9,
#     max_target_seqs = 500,
#     out_dir = "exports/blastn_pass2"
#   )
# 
#   res_pass2$db_source <- "full_reference"
# } else {
#   res_pass2 <- res_pass1[0, ]
# }
# 
# rb_disconnect(con)
# 
# # combine both assignment tabls
# final_assignments <- bind_rows(
#   res_pass1 |>
#     filter(assignment_status != "no_match") |>
#     mutate(db_source = "12S_subset"),
#   res_pass2
# )
# 
# missing_ids <- setdiff(all_asv_ids, final_assignments$asv_id)
# 
# if (length(missing_ids) > 0) {
#   pad <- final_assignments[rep(NA_integer_, length(missing_ids)), , drop = FALSE]
#   pad$asv_id <- missing_ids
#   pad$assignment_status <- "no_match"
#   pad$db_source <- "no_hit"
#   final_assignments <- bind_rows(final_assignments, pad)
# }
# 
# head(final_assignments)
# 
# # export taxonomic assignment table
# write_csv(
#   final_assignments,
#   "exports/sava_river_12s_assignments_two_pass.csv"
# )

## ----eval=TRUE----------------------------------------------------------------
head(final_assignments)

## -----------------------------------------------------------------------------
# community_matrix <- rb_build_species_matrix(
#   assignment = final_assignments,
#   asv_table = sava_counts,
#   include_ambiguous = FALSE,
#   unassigned = FALSE
# )
# 
# community_matrix[1:5, 1:5]
# 
# # export species-site matrix
# write_csv(
#   community_matrix,
#   "exports/sava_river_community_matrix.csv"
# )

## ----eval=TRUE----------------------------------------------------------------
community_matrix[1:5, 1:5]

## ----eval=TRUE----------------------------------------------------------------
PrintBibliography(myBib)

