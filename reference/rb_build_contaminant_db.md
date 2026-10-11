# Build contaminant reference databases for QC screening

Functions to download contaminant sequences from NCBI and build
nucleotide BLAST databases for screening reference sequences against
known contaminants (e.g., human mtDNA, lab bacteria, NUMTs).

- `rb_build_contaminant_db()`: Downloads sequences from NCBI based on an
  organisms table, saving one FASTA file per category.

- `rb_build_blast_db()`: Combines FASTA files and builds a nucleotide
  BLAST database using `makeblastdb`. Creates an intermediate combined
  FASTA file (e.g., `prefix_combined.fasta`) in the same directory as
  the `out_prefix`. This file is retained on disk after BLAST database
  is built.

## Usage

``` r
rb_build_contaminant_db(
  organisms,
  out_dir,
  retmax = 200,
  requests_per_second = 3
)

rb_build_blast_db(
  fasta_paths,
  out_prefix,
  title = "Contaminant Database",
  makeblastdb = "makeblastdb",
  stream_combine = FALSE
)
```

## Arguments

- organisms:

  A data frame with columns:

  - `type`: Either `"accession"` (direct `entrez_fetch`) or `"query"`
    (`entrez_search` then fetch).

  - `category`: Output grouping that becomes the FASTA filename (e.g.,
    `"human_mtDNA"`, `"lab_bacteria"`, `"numt_source"`).

  - `value`: The accession ID (for `type = "accession"`) or NCBI search
    query string (for `type = "query"`).

- out_dir:

  Directory where output FASTA files will be saved. Created
  automatically if it does not exist.

- retmax:

  Maximum number of sequences to retrieve per NCBI search.

- requests_per_second:

  Throttle rate between NCBI requests. The unauthenticated E-utilities
  limit is approximately 3 requests per second.

- fasta_paths:

  Character vector of FASTA file paths to combine and build into a BLAST
  database.

- out_prefix:

  Output prefix for the BLAST database files (e.g.,
  `"data/contaminants/contaminant_db"`).

- title:

  Title string embedded in the BLAST database metadata.

- makeblastdb:

  Path to the `makeblastdb` executable.

- stream_combine:

  Logical. If `TRUE`, combines files via chunked streaming instead of
  loading them fully into memory. Recommended for large contaminant
  sets.

## Value

- `rb_build_contaminant_db()` returns a named character vector of output
  FASTA file paths, one per category.

- `rb_build_blast_db()` returns the `out_prefix` path invisibly.

## Details

The typical workflow is:

1.  Define an `organisms` table specifying what to download.

2.  Call `rb_build_contaminant_db()` to download FASTA files.

3.  Call `rb_build_blast_db()` on the relevant FASTA files to create the
    BLAST database.

Note: NUMT sequences are typically screened via exact substring matching
([`rb_screen_numts()`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md))
rather than BLAST, so NUMT FASTA files are usually excluded from the
BLAST database.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define organisms to download
organisms <- data.frame(
  type = c("accession", "query"),
  category = c("human_mtDNA", "lab_bacteria"),
  value = c(
    "NC_012920.1",
    "(Escherichia coli[ORGN] OR Pseudomonas[ORGN]) AND 16S[GENE]"
  ),
  stringsAsFactors = FALSE
)

# Download contaminant sequences
fasta_paths <- rb_build_contaminant_db(organisms, out_dir = "data/contaminants")

# Build BLAST database (excluding NUMTs if present)
blast_paths <- fasta_paths[!grepl("numt", names(fasta_paths))]
rb_build_blast_db(blast_paths, out_prefix = "data/contaminants/contaminant_db")
} # }
```
