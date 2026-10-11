# Curating an eDNA Reference Database for Slovenian Freshwater Fishes

## Introduction

`barcurateR` provides a reproducible workflow for building, auditing,
and exporting curated DNA barcode reference databases for eDNA
metabarcoding, which includes basic taxonomic assignment functions,
i.e., BLASTn and exact matching.

This vignette walks through a realistic freshwater fish case study:
curating a local reference database using NCBI and BOLD sequences to be
used for the taxonomic assignment for a eDNA dataset.

### Dataset

For this example, we used the Slovenian fish species list presented by
the Slovenian Native Fish Society based on the “Freshwater fish and
lampreys in Slovenia” (Povž et al., 2015; Društvo za preučevanje rib
Slovenije, 2018), supplemented by “The Updated List of Non-Native
Freshwater Fishes in Slovenia with Note of their Potential Impacts in
Inland Waters” (Povž et al., 2018), for a total of 97 species. The
sequences from NCBI (queried and retrieved on 2026-10-02) and BOLD
(retrieved on 2026-09-16) were parsed and subsetted separately based on
the compiled species list.

The reference database validation was conducted on fish eDNA dataset
from the middle reach of the Sava River as reported from the study:
“Temporal and Spatial Dynamics of the Lotic Fish Communities: A
Comparison of Coffee Filter-Based Passive eDNA Collection Versus Active
eDNA Filtering” (Mlinarec et al., 2025). The validation eDNA ASV was
reconstructed following the DADA2 workflow (Callahan et al., 2016) using
parameters described by Mlinarec et al. (2025), where applicable.

*Note: File names and paths may need to be replaced when processing your
own data; curation and taxonomic assignment steps may require external
software i.e., BLAST, MAFFT, and FastTree.*

## 1. Setup and Configuration

### Loading libraries

`dplyr`, `reader`, and `tidyr` are used for data manipulation and
`Biostrings` for FASTA handling.

``` r

library(barcurateR)
library(dplyr)
library(readr)
library(tidyr)
library(Biostrings)
```

### Loading the Dataset

For this case study, reference sequences retrieved from NCBI
(`ncbi_raw`) and BOLD (`bold_raw`), along with a taxonomic table
(`taxonomy_table`) for the target species list.

``` r

ncbi_csv <- "data-raw/slovenian_fish/raw/ncbi_sequence_raw.csv"
bold_csv <- "data-raw/slovenian_fish/raw/bold_sequence_raw.csv"
taxonomy_csv <- "data-raw/slovenian_fish/processed/slovenian_fish_taxonomy.csv"

ncbi_raw <- read_csv(ncbi_csv, show_col_types = FALSE) |>
  filter(!is.na(sequence))

bold_raw <- read_csv(bold_csv, show_col_types = FALSE) |>
  filter(!is.na(sequence))

taxonomy_table <- read_csv(taxonomy_csv, show_col_types = FALSE)
```

## 2. Standardization & ambiguity resolution

### Data standardization

We use
[`rb_parse_source_table()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md)
to standardize the heterogeneous inputs. This function extracts sequence
marker information from the sequence headers (common in NCBI/BOLD
downloads).
[`rb_parse_source_table()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md)
allows users to provide the column mapping and extracts additional
marker information from the `description_col` as typically included in
the header from NCBI and BOLD data. Detailed descriptions for required
variables are described in the manual for
[`rb_parse_source_table()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md).

``` r

# NCBI
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

# BOLD
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
```

### Data merge

`rb_combine_source()` automatically merges standardized sequence data
and deduplicates identical species-sequence entries.

``` r

combined <- rb_combine_sources(list(ncbi_parsed, bold_parsed))
```

### Ambiguity resolution

Often, sequences downloaded from large public repositories, e.g., NCBI,
also contain identical sequences attributed to different taxa.
[`rb_detect_ambiguity()`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md)
helps users identify these sequences by flagging many-to-one
species-sequence data for subsequent action and exports the ambiguities
as a CSV as documentation. Users can provide an ambiguity resolution
table for
[`rb_resolve_ambiguous()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md),
confirm that the sequence dataset no longer has any ambiguities
(`"stop"` or `"warn"`), or as we show here, to automatically drop all
ambiguous sequences.

``` r

combined <- rb_detect_ambiguity(combined)
n_ambig <- sum(combined$is_ambiguous)

if (n_ambig > 0) {
  rb_write_ambiguity_csv(
    combined[combined$is_ambiguous, ],
    "data-raw/slovenian_fish/processed/ambiguous_sequences.csv"
  )
}

resolved <- rb_resolve_ambiguous(
  combined,
  on_unresolved = "drop"
)
```

## 3. Curation Pipeline

The curation workflow for `barcurateR` reproduces the pipeline described
from “YZFishDB: A curated DNA reference database for the fishes of the
Yangtze River Basin” . the core function
[`rb_curate_reference()`](https://zkua01.github.io/barcurateR/reference/rb_curate_reference.md)
completes the entire curation workflow as described by Kua et al. (2026)
and exports as SQLite database, while allowing users to opt for major
steps in the workflow. Specific requirements and dependencies for
[`rb_curate_reference()`](https://zkua01.github.io/barcurateR/reference/rb_curate_reference.md)
are detailed in the manual (`?rb_curate_reference()`).

For this case study, we passed our combined sequence dataset through the
full curation workflow by providing the contaminant `BLASTn` database
and the numts `FASTA`.

*Note: External software BLAST+ (`blastn` and `makeblastdb`), MAFFT, and
FastTree are required to complete the full curation pipeline*

``` r

contam_path <- "tests/testthat/fixtures/contaminants/contaminant_db"
numts_path  <- "tests/testthat/fixtures/contaminants/fish_numts.fasta"
db_path <- "data-raw/slovenian_fish/processed/slovenian_fish_refdb.sqlite"

curated <- rb_curate_reference(
  data = resolved,
  blast_db = contam_path,
  numt_fasta = numts_path,
  run_classifier = TRUE,
  run_divergence = TRUE,  # requires MAFFT and FastTree
  run_barcode_gap = TRUE,
  taxonomy_table = taxonomy_table,
  db_path = db_path,
  on_ambiguous = "stop",
  require_taxonomy = TRUE,
  archive_ambiguity = TRUE,
  ambiguity_path = "data-raw/slovenian_fish/processed/ambiguous_sequences.csv"
)
```

## 4. Database auditing

The curated reference database can be queried and audited through
several built-in functions to inspect and summarize information e.g.,
taxonomic and marker coverage. These functions allow users to better
understand the composition, quality, and coverage of the reference
database. Curated `barcurateR` reference database can be connected
through
[`rb_connect()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md).

``` r

con <- rb_connect(db_path)
```

### Reference database sequence quality and coverage

#### Audit QC flags

``` r

qc_sum <- rb_qc_summary(con)
print(qc_sum)
```

``` r

print(outputs[[1]])
#>                                           qc_flag n_sequences
#> 1                                            pass       10431
#> 2                            manual_review_needed        5902
#> 3                                    frameshifted        3317
#> 4 tree_high_length_variation|manual_review_needed        1850
#> 5                  divergent|manual_review_needed         332
#> 6                                   internal_stop         284
#> 7                      internal_stop|frameshifted          80
#> 8                                    sequence_gap          31
#> 9                                  sequence_short          13
```

#### Audit marker coverage

``` r

marker_cov <- rb_marker_coverage(con, rank = "order")
head(marker_cov)
```

``` r

head(outputs[[2]])
#>              order seq_type n_sequences n_sources
#> 1 Acipenseriformes      12S          11         1
#> 2 Acipenseriformes      16S           7         1
#> 3 Acipenseriformes      18S           1         1
#> 4 Acipenseriformes      COI          60         2
#> 5 Acipenseriformes     CYTB          68         2
#> 6 Acipenseriformes   genome          10         1
```

#### Audit taxonomic coverage

``` r

taxa_cov <- rb_list_taxa(con, rank = "species")
head(taxa_cov)
```

``` r

head(outputs[[3]])
#>                   species n_sequences
#> 1           Abramis brama          60
#> 2        Acipenser baerii          69
#> 3      Acipenser ruthenus          53
#> 4 Alburnoides bipunctatus         132
#> 5       Alburnus alburnus          77
#> 6      Alburnus arborella          35
```

### Barcode gap analysis

To help users evaluate the efficacy of taxonomic assignments,
`barcurateR` includes a simple barcode gap analysis.

``` r

gap_metrics <- rb_barcode_gap(con)

if (nrow(gap_metrics) > 0) {
  poor_resolution <- gap_metrics |>
    filter(taxonomic_resolution != "species" | assignment_risk == "high")
}

head(poor_resolution)
```

``` r

head(outputs[[4]])
#>                   species marker n_sequences max_intra    min_inter gap_exists
#> 1           Abramis brama   CYTB          23 0.3126521 0.0030413625          0
#> 2        Acipenser baerii   CYTB          25 0.3111314 0.0009124088          0
#> 3      Acipenser ruthenus   CYTB          27 0.3254258 0.0072992701          0
#> 4 Alburnoides bipunctatus   CYTB         105 0.0586983 0.0413625304          0
#> 5       Alburnus alburnus   CYTB          24 0.3086983 0.0036496350          0
#> 6      Alburnus arborella   CYTB          31 0.1806569 0.0082116788          0
#>     gap_width p_overlap q_overlap median_intra median_inter   intra_95
#> 1 -0.30961071    0.7549    0.8848  0.105535280    0.1517640 0.30748175
#> 2 -0.31021898    0.9500    0.9467  0.174574209    0.2010341 0.26434002
#> 3 -0.31812652    0.7550    0.9436  0.137773723    0.1563260 0.31219586
#> 4 -0.01733577    0.4762    0.0235  0.012165450    0.1146594 0.05687348
#> 5 -0.30504866    0.8043    0.9535  0.046836983    0.1322993 0.28094586
#> 6 -0.17244526    0.4860    0.7144  0.007603406    0.1055353 0.16545012
#>      inter_5 percentile_gap robust_gap n_congeners_in_db min_congeneric
#> 1 0.06082725   -0.246654501          0                 0             NA
#> 2 0.07968370   -0.184656326          0                 1   0.0009124088
#> 3 0.09124088   -0.220954988          0                 1   0.0009124088
#> 4 0.06356448    0.006690998          1                 0             NA
#> 5 0.06569343   -0.215252433          0                 1   0.0082116788
#> 6 0.06265207   -0.102798054          0                 1   0.0082116788
#>   congeneric_gap congeneric_p_overlap congeneric_q_overlap n_diagnostic_sites
#> 1             NA                   NA                   NA                  0
#> 2     -0.3102190               0.9500               0.9926                  0
#> 3     -0.3245134               0.8205               0.9985                  0
#> 4             NA                   NA                   NA                  0
#> 5     -0.3004866               0.6884               0.9906                  0
#> 6     -0.1724453               0.4860               0.7782                  0
#>   sampling_adequacy taxonomic_resolution gap_may_be_undersampling_artifact
#> 1               low      genus_or_higher                                 0
#> 2          moderate      genus_or_higher                                 0
#> 3          moderate      genus_or_higher                                 0
#> 4               low      genus_or_higher                                 0
#> 5          moderate      genus_or_higher                                 0
#> 6          moderate      genus_or_higher                                 0
#>   error_msg assignment_risk recommended_threshold
#> 1      <NA>            high                  0.97
#> 2      <NA>            high                  0.97
#> 3      <NA>            high                  0.97
#> 4      <NA>            high                  0.97
#> 5      <NA>            high                  0.97
#> 6      <NA>            high                  0.97
```

``` r

rb_disconnect(con)
```

## 5. Downstream Applications

### Exporting database

For users conducting most of the bioinformatics and downstream taxonomic
assignment work in other popular platforms, `barcurateR` provides a
simple wrapper function to export the selected reference database for i)
DADA2, ii) QIIME2, iii) BLASTn, iv) USEARCH/VSEARCH. Here we first
subset the curated database based on the `marker` and `qc_flag`
according to our needs using
[`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
before exporting them.

``` r

con <- rb_connect(db_path)

refs_12S <- rb_get_sequences(con, marker = "12S", qc_flag = "pass")
refs_COI <- rb_get_sequences(con, marker = "COI", qc_flag = "pass")

# DADA2
rb_export_dada2(
  refs_12S,
  "exports/slovenian_fish_12S_dada2.fasta"
)

# QIIME2
rb_export_qiime2(
  refs_COI,
  "exports/slovenian_fish_COI_qiime2.fasta",
  "exports/slovenian_fish_COI_taxonomy.tsv"
)

# BLASTn: 12S only
rb_export_blastn(
  rb_get_sequences(con, marker = "12S", qc_flag = "pass"),
  "exports/slovenian_fish_12s_blastn.fasta"
)

# BLASTn: all markers
rb_export_blastn(
  rb_get_sequences(con, qc_flag = "pass"),
  "exports/slovenian_fish_all_markers_blastn.fasta"
)

rb_disconnect(con)
```

### Taxonomic assignment

`barcurateR` includes two taxonomic assignment methods, i.e., BLASTn and
exact matching. BLASTn taxonomic assignment relies on the external
software from the BLAST+ suite (`blastn`) while exact matching is a
dependency-free method that was included for testing and quick
assignment jobs.

The package requires ASV sequences in a FASTA file (`sava_asvs_fasta`)
and a corresponding count table (`sava_counts`).

``` r

sava_asvs_fasta <- "data-raw/sava_river/processed/sava_asvs.fasta"
sava_counts <- read_csv(
  "data-raw/sava_river/processed/sava_abundance_table.csv",
  show_col_types = FALSE
)
```

For this example, we conducted a two-pass taxonomic assignment using
BLASTn:

Pass 1: Assign against the 12S subset of the curated reference database
(`marker="12S"` and `qc_flag="pass"`). Pass 2: Assign unassigned ASVs
against the full reference database (`marker=NULL` and
`qc_flag=c("pass","manual_review_needed")`) to rescue borderline
matches.

``` r

con <- rb_connect(db_path)

# Pass 1: 12S subset
res_pass1 <- rb_assign_edna(
  asv_fasta = sava_asvs_fasta,
  con = con,
  marker = "12S",
  method = "blastn",
  min_identity = 99,
  out_dir = "exports/blastn_pass1"
)

all_asv_ids <- names(Biostrings::readDNAStringSet(sava_asvs_fasta))
assigned_ids <- res_pass1$asv_id[res_pass1$assignment_status != "no_match"]
unassigned_ids <- setdiff(all_asv_ids, assigned_ids)

# Pass 2: full reference for unassigned ASVs
if (length(unassigned_ids) > 0) {
  all_seqs <- Biostrings::readDNAStringSet(sava_asvs_fasta)
  pass2_fasta <- tempfile(fileext = ".fasta")

  Biostrings::writeXStringSet(
    all_seqs[names(all_seqs) %in% unassigned_ids],
    pass2_fasta
  )

  res_pass2 <- rb_assign_edna(
    asv_fasta = pass2_fasta,
    con = con,
    marker = NULL,
    method = "blastn",
    qc_flag = c("pass", "manual_review_needed"),
    min_identity = 99,
    min_coverage = 0.9,
    max_target_seqs = 500,
    out_dir = "exports/blastn_pass2"
  )

  res_pass2$db_source <- "full_reference"
} else {
  res_pass2 <- res_pass1[0, ]
}

rb_disconnect(con)

# combine both assignment tabls
final_assignments <- bind_rows(
  res_pass1 |>
    filter(assignment_status != "no_match") |>
    mutate(db_source = "12S_subset"),
  res_pass2
)

missing_ids <- setdiff(all_asv_ids, final_assignments$asv_id)

if (length(missing_ids) > 0) {
  pad <- final_assignments[rep(NA_integer_, length(missing_ids)), , drop = FALSE]
  pad$asv_id <- missing_ids
  pad$assignment_status <- "no_match"
  pad$db_source <- "no_hit"
  final_assignments <- bind_rows(final_assignments, pad)
}

head(final_assignments)

# export taxonomic assignment table
write_csv(
  final_assignments,
  "exports/sava_river_12s_assignments_two_pass.csv"
)
```

``` r

head(final_assignments)
#> # A tibble: 6 × 16
#>   asv_id species         genus family pident  qcov scovs ccovs   evalue bitscore
#>   <chr>  <chr>           <chr> <chr>   <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 ASV_3  Squalius cepha… Squa… Cypri…  100   1     0.931 0.965 1.39e-84      302
#> 2 ASV_5  Alburnus albur… Albu… Cypri…  100   1     0.931 0.965 1.39e-84      302
#> 3 ASV_7  Barbus barbus   Barb… Cypri…  100   1     0.931 0.965 4.98e-84      300
#> 4 ASV_9  Vimba vimba     Vimba Cypri…  100   1     0.931 0.965 1.39e-84      302
#> 5 ASV_12 Abramis brama   Abra… Cypri…  100   1     0.932 0.965 3.90e-85      303
#> 6 ASV_13 Chondrostoma s… Chon… Cypri…   99.4 0.994 0.926 0.959 2.33e-82      294
#> # ℹ 6 more variables: bitscore_pb <dbl>, n_best_species <dbl>,
#> #   assignment_status <chr>, confidence <chr>, tie_taxa <lgl>, db_source <chr>
```

After the taxonomic assignment is completed, we can use
[`rb_build_species_matrix()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md)
to create the species-site data table using the taxonomic assignment
(`final_assignments`) and the count table (`sava_counts`).

``` r

community_matrix <- rb_build_species_matrix(
  assignment = final_assignments,
  asv_table = sava_counts,
  include_ambiguous = FALSE,
  unassigned = FALSE
)

community_matrix[1:5, 1:5]

# export species-site matrix
write_csv(
  community_matrix,
  "exports/sava_river_community_matrix.csv"
)
```

``` r

community_matrix[1:5, 1:5]
#> # A tibble: 5 × 5
#>   species             SRR31349521 SRR31349522 SRR31349523 SRR31349524
#>   <chr>                     <dbl>       <dbl>       <dbl>       <dbl>
#> 1 Abramis brama               251         187          38          42
#> 2 Alburnus alburnus           511         272          59          55
#> 3 Ballerus ballerus          1203        1013          39          30
#> 4 Barbatula barbatula           0           0           0           0
#> 5 Barbus barbus               764         610         146         118
```

## Conclusion

This case study demonstrates how `barcurateR` streamlines the complex
process of reference database curation. By packaging robust QC,
ambiguity resolution, and auditing tools, `barcurateR` ensures that
downstream eDNA assignments are based on high-quality, transparent, and
reproducible reference data.

## References

\[1\] B. J. Callahan et al. “DADA2: High-resolution sample inference
from Illumina amplicon data”. En. In: *Nature Methods* 13.7 (Jul. 2016),
pp. 581–583. ISSN: 1548-7091, 1548-7105. DOI:
[10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869). URL:
<https://www.nature.com/articles/nmeth.3869> (visited on 08/07/2025).

\[2\] Društvo za preučevanje rib Slovenije. *Naše ribe*. Si. Ljubljana,
2018. URL: <https://www.dprs.si/sl/ribe-v-sloveniji> (visited on
09/16/2026).

\[3\] Z. X. Kua et al. “YZFishDB: A curated DNA reference database for
the fishes of the Yangtze River Basin”. En. In: *Aquatic Diversity and
Ecology* 1.1 (Jun. 2026), p. 100009. ISSN: 31178634. DOI:
[10.1016/j.ade.2026.100009](https://doi.org/10.1016/j.ade.2026.100009).
URL: <https://linkinghub.elsevier.com/retrieve/pii/S3117863426000092>
(visited on 06/13/2026).

\[4\] J. Mlinarec et al. “Temporal and Spatial Dynamics of the Lotic
Fish Communities: A Comparison of Coffee Filter‐Based Passive eDNA
Collection Versus Active eDNA Filtering”. En. In: *Environmental DNA*
7.2 (Mar. 2025), p. e70065. ISSN: 2637-4943, 2637-4943. DOI:
[10.1002/edn3.70065](https://doi.org/10.1002/edn3.70065). URL:
<https://onlinelibrary.wiley.com/doi/10.1002/edn3.70065> (visited on
09/15/2026).

\[5\] M. Povž et al. *Freshwater fish and lampreys in Slovenia \[in
Slovene\]*. Si. Ljubljana: Zavod Umbra, 2015.

\[6\] M. Povž et al. “The Updated List of the Non-Native Freshwater
Fishes in Slovenia with Note of their Potential Impact in Inland
Waters”. En. In: *Pakistan Journal of Zoology Supplementary Series* 13
(2018), pp. 1–7.
