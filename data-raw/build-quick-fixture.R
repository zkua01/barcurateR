# ============================================================
# data-raw/build_quick_fixture.R (Robust Version)
# ============================================================

# 1. Dependencies
pkgs <- c("rentrez", "bold", "readr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(rentrez)
library(bold)
library(readr)

# 2. Configuration
species_list <- c(
  "Alligator sinensis",
  "Python molurus",
  "Umbra krameri",
  "Esox lucius"
)

max_per_marker <- 2
fixture_dir <- file.path("tests", "testthat", "fixtures")
contaminant_dir <- file.path(fixture_dir, "contaminants")
dir.create(contaminant_dir, recursive = TRUE, showWarnings = FALSE)

# Optional NCBI API key
if (nzchar(Sys.getenv("ENTREZ_API_KEY"))) {
  rentrez::set_entrez_key(Sys.getenv("ENTREZ_API_KEY"))
  delay <- 0.15
} else {
  delay <- 0.4
}

# 3. NCBI Download
ncbi_terms <- list(
  "12S" = '(12S[Title] OR "12S ribosomal RNA"[Title] OR "MT-RNR1"[Gene Name])',
  "16S" = '(16S[Title] OR "16S ribosomal RNA"[Title] OR "MT-RNR2"[Gene Name])',
  "18S" = '(18S[Title] OR "18S ribosomal RNA"[Title])',
  "COI" = '(COI[Title] OR CO1[Title] OR "cytochrome c oxidase"[Title] OR "MT-CO1"[Gene Name])' # Fallback if BOLD fails
)

parse_fasta <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  lines <- strsplit(txt, "\r?\n")[[1]]
  headers <- sub("^>", "", lines[grepl("^>", lines)])
  seqs <- character(length(headers))
  idx <- grep("^>", lines)
  for (i in seq_along(idx)) {
    start <- idx[i] + 1
    end <- ifelse(i == length(idx), length(lines), idx[i+1] - 1)
    if (end >= start) seqs[i] <- paste(lines[start:end], collapse = "")
  }
  data.frame(
    sequence_id = sub("\\s+.*", "", headers),
    description = headers,
    sequence = toupper(gsub("\\s+", "", seqs)),
    stringsAsFactors = FALSE
  )
}

ncbi_results <- list()
for (sp in species_list) {
  for (mk in names(ncbi_terms)) {
    cat(sprintf("NCBI: %s / %s\n", sp, mk))
    term <- paste0(sp, "[Organism] AND ", ncbi_terms[[mk]])
    
    search <- tryCatch(
      rentrez::entrez_search(db = "nuccore", term = term, retmax = max_per_marker),
      error = function(e) list(ids = character(0))
    )
    
    if (length(search$ids) > 0) {
      Sys.sleep(delay)
      fasta <- tryCatch(
        rentrez::entrez_fetch(db = "nuccore", id = search$ids, rettype = "fasta"),
        error = function(e) NULL
      )
      parsed <- parse_fasta(fasta)
      if (!is.null(parsed)) {
        parsed$source <- "ncbi"
        parsed$species_query <- sp
        parsed$marker_query <- mk
        ncbi_results[[length(ncbi_results) + 1]] <- parsed
      }
    } else {
      cat("  -> No results\n")
    }
    Sys.sleep(delay)
  }
}

ncbi_tbl <- do.call(rbind, ncbi_results)
if (is.null(ncbi_tbl)) ncbi_tbl <- data.frame()
write_csv(ncbi_tbl, file.path(fixture_dir, "ncbi_marker_sequences.csv"))
cat(sprintf("\nSaved %d NCBI sequences.\n", nrow(ncbi_tbl)))


# 4. BOLD data from local source

bold_data <- data.frame(
  processid = c("ANGBF35881-19", "ANGBF5754-12"),
  sampleid = c("JX261985", "HM563702"),
  fieldid = c("", "HM563702"),
  museumid = c("", ""),
  record_id = c("ANGBF35881-19.COI-5P", "ANGBF5754-12.COI-5P"),
  specimenid = c("10134239", "3017389"),
  bin_uri = c("BOLD:AAO6269", "BOLD:AAA5988"),
  taxid = c("342560", "29373"),
  kingdom = c("Animalia", "Animalia"),
  phylum = c("Chordata", "Chordata"),
  class = c("Actinopterygii", "Actinopterygii"),
  order = c("Esociformes", "Esociformes"),
  family = c("Umbridae", "Esocidae"),
  genus = c("Umbra", "Esox"),
  species = c("Umbra krameri", "Esox lucius"),
  identification_rank = c("species", "species"),
  marker_code = c("COI-5P", "COI-5P"),
  nuc = c(
    paste0(
      "AAAGACATTGGCACCCTTTATTTAGTATTTGGTGCCTGGGCCGGAATAGTCGGCACCGCTCTAAGCCTTTTGATTCGGGCTGAACTT",
      "AGCCAGCCGGGGGCCTTGCTCGGAGATGACCAAATTTATAATGTTATCGTCACTGCACACGCCTTTGTCATAATTTTCTTTATAGTA",
      "ATGCCCATTATAATTGGAGGCTTTGGAAACTGATTAGTCCCCCTTATGATTGGGGCCCCAGACATAGCATTTCCTCGAATAAATAAT",
      "ATAAGCTTCTGGCTCCTTCCTCCTTCTTTCCTTCTCCTCCTAGCCTCTTCAGGGGTTGAAGCAGGGGCAGGAACAGGATGAACTGTC",
      "TATCCCCCTCTCGCCGGTAACCTTGCCCATGCAGGCGCTTCCGTAGACTTAACTATTTTTTCTCTACATTTGGCTGGAGTTTCCTCA",
      "ATTCTGGGCGCCATCAATTTTATTACAACAATTATCAACATAAAACCCCCTGCAATCTCACAATACCAAACACCCCTATTTGTCTGG",
      "GCAGTCCTAATCACAGCAGTCCTTTTATTACTTTCTCTTCCCGTTCTCGCCGCAGGGATCACAATGCTACTCACGGATCGAAACCTA",
      "AACACCACCTTTTTTGACCCAGCCGGAGGAGGAGACCCTATTTTATATCAACACCTTTTTTGATTCTTTGGCCACCC"
    ),
    paste0(
      "CTTTATTTAGTATTTGGTGCTTGAGCCGGAATAGTCGGCACAGCCTTAAGCCTTTTAATCCGGGCCGAACTAAGCCAGCCAGGGGCT",
      "CTCTTAGGTGACGACCAGATTTATAATGTTATCGTTACAGCCCATGCCTTTGTTATAATCTTTTTTATAGTTATACCCGTTATAATT",
      "GGGGGTTTTGGAAACTGATTAATTCCCCTAATGATTGGTGCCCCCGACATGGCCTTCCCCCGCATAAATAATATAAGCTTCTGACTT",
      "CTCCCCCCCTCCTTTTTACTTCTCTTGGCCTCCTCAGGTGTTGAAGCTGGTGCTGGTACTGGCTGAACAGTTTATCCGCCTTTGGCC",
      "GGAAACTTAGCACACGCAGGTGCTTCTGTAGACTTAACTATTTTCTCTCTCCACCTGGCCGGAATTTCTTCTATTCTAGGAGCTATT",
      "AATTTTATTACCACAATTATTAACATAAAACCCCCCGCCATCTCACAATATCAGACACCATTATTTGTTTGAGCAGTCCTGATTACA",
      "GCTGTACTTCTACTTCTATCTCTCCCAGTCCTAGCCGCTGGAATTACCATATTGCTCACAGACCGAAATTTAAACACCACATTCTTT",
      "GACCCCGCTGGTGGTGGAGACCCTATTCTATACCAACACCTC"
    )
  ),
  nuc_basecount = c(686, 651),
  insdc_acs = c("JX261985", "HM563702"),
  stringsAsFactors = FALSE
)

readr::write_tsv(
  bold_data,
  "tests/testthat/fixtures/bold_coi_sequences.tsv"
)

message("Created BOLD fixture with ", nrow(bold_data), " sequences.")

# 5. Contaminants (Human mtDNA)
cat("\nDownloading Human mtDNA contaminant...\n")
human_fasta <- tryCatch(
  rentrez::entrez_fetch(db = "nuccore", id = "NC_012920.1", rettype = "fasta"),
  error = function(e) NULL
)
if (!is.null(human_fasta)) {
  writeLines(human_fasta, file.path(contaminant_dir, "human_mtDNA.fasta"))
  cat("Saved human_mtDNA.fasta\n")
}

# 6. Taxonomy Fixture
tax_tbl <- data.frame(
  species = species_list,
  genus = c("Alligator", "Python", "Umbra", "Esox"),
  family = c("Alligatoridae", "Pythonidae", "Umbridae", "Esocidae"),
  order = c("Crocodilia", "Squamata", "Salmoniformes", "Salmoniformes"),
  class = c("Reptilia", "Reptilia", "Actinopterygii", "Actinopterygii"),
  phylum = rep("Chordata", 4),
  kingdom = rep("Animalia", 4),
  stringsAsFactors = FALSE
)
write_csv(tax_tbl, file.path(fixture_dir, "quick_taxonomy.csv"))

cat("\n=== FIXTURE BUILD COMPLETE ===\n")