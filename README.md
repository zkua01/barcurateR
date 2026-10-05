# barcurateR

`barcurateR` is an R package designed to **standardize, curate, audit, and deploy DNA reference databases** for ecological metabarcoding. It combines the data curation workflow developed for YZFishDB with the database auditing and taxonomic assignment functions from `regionbarcoder`, providing a complete, reproducible pipeline from raw sequence data to taxonomic assignment.

## Installation

You can install the development version of `barcurateR` from GitHub using the `remotes` or `pak` package:

```{r}
# install.packages("remotes")
remotes::install_github("zkua01/barcurateR")
```

# The `barcurateR` Workflow

The package is built around a modular, 5-step workflow.

## 1. Parse and Standardize Sources

Import raw sequence tables from various sources (e.g., NCBI, BOLD) and standardize them into a unified schema. If a description column is provided, barcurateR will automatically infer marker types (e.g., COI, 12S) and completeness.

```{r}
parsed_ncbi <- rb_parse_source_table(
  raw_data,
  source_name = "ncbi",
  column_map = c(species = "scientific_name", sequence = "nt_seq"),
  description_col = "definition"
)
```

## 2. Combine and Resolve Ambiguities

Merge standardized tables from multiple sources and identify sequences assigned to multiple species. Ambiguities must be resolved before Quality Control (QC) to ensure reliable downstream metrics.

```{r}
combined <- rb_combine_sources(list(ncbi_parsed, bold_parsed))
combined <- rb_detect_ambiguity(combined)
```

## 3. Run the Curation Pipeline

Execute the comprehensive QC wrapper. This step screens for contaminants and NUMTs, checks codon/rRNA integrity, optionally runs ML-based marker classification, joins taxonomic lineages, and writes the curated data to a SQLite database.

```{r}
curated <- rb_curate_reference(
  data = resolved,
  blast_db = "path/to/contaminant_db",
  numt_fasta = "path/to/numts.fasta",
  taxonomy_table = tax_lineage,
  run_classifier = TRUE,
  db_path = "curated_reference.sqlite"
)
```

## 4. Query, Audit, and Export

Interact with your curated SQLite database to audit data quality, check marker coverage, and export sequences in formats compatible with popular metabarcoding pipelines (BLAST, DADA2, QIIME2, USEARCH).

```{r}
con <- rb_connect("curated_reference.sqlite")

# Audit the database
rb_qc_summary(con)
rb_marker_coverage(con, rank = "genus")

# Export for DADA2
refs <- rb_get_sequences(con, qc_flag = "pass")
rb_export_dada2(refs, "reference_dada2.fasta")

rb_disconnect(con)
```

## 5. Assign eDNA Sequences

Assign Amplicon Sequence Variants (ASVs) to your curated reference database using either dependency-free exact matching or BLASTn.

```{r}
assignmentt <- rb_assign_edna(
  asv_fasta = "asvs.fasta",
  db_path = "curated_reference.sqlite",
  marker = "12S",
  method = "blastn", # or method = "exact"
  min_identity = 99
)
```

# External Software Requirements

While core parsing, QC, and exact-matching functions run natively in R, several advanced features require external command-line tools to be installed and available on your system PATH:

  - BLAST+ (`blastn`, `makeblastdb`): Required for contaminant screening (`rb_screen_contaminants` or `rb_curate_reference(blast_db = "path/to/contaminant.db")`) and BLAST-based taxonomic assignment (`method = "blastn"`).
  
  - MAFFT and FastTree: Required for phylogenetic divergence checks (`run_divergence = TRUE` in `rb_curate_reference()`).
  
A more detailed example workflow using for assigning a Sava River ASV dataset is available as a vignette (`browseVignettes("barcurateR")`).

For questions or bug reports, please use the [GitHub Issues page](https://github.com/zkua01/barcurateR/issues).