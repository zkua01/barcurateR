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

