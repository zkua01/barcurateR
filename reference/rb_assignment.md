# Assign eDNA sequences to reference taxa

Functions to assign amplicon sequence variants (ASVs) to a curated
reference database using either exact sequence matching or BLASTn
alignment.

- `rb_assign_edna()`: Main wrapper for taxonomic assignment.

- `rb_assign_blastn()`: Runs BLASTn against a reference database.

- `rb_exact_hits()`: Finds exact sequence matches in the reference.

- `rb_score_assignments()`: Resolves and scores assignment hits,
  handling ties.

- `rb_parse_blast_tabular()`: Parses standard BLAST tabular output.

- `rb_build_species_matrix()`: Aggregates ASV assignments into a
  species-by-sample matrix.

- `rb_read_fasta()`: Reads a FASTA file into a data frame (used
  internally for ASV input).

## Usage

``` r
rb_read_fasta(path)

rb_score_assignments(
  hits,
  min_identity = 99,
  min_coverage = 0.9,
  tie_bitscore = 0
)

rb_exact_hits(asvs, refs)

rb_parse_blast_tabular(path)

rb_assign_blastn(
  asv_fasta,
  refs,
  out_dir,
  min_identity,
  min_coverage,
  max_target_seqs = 100,
  blastn,
  makeblastdb
)

rb_assign_edna(
  asv_fasta,
  con = NULL,
  db_path = NULL,
  marker = NULL,
  qc_flag = "pass",
  method = c("blastn", "exact"),
  min_identity = 99,
  min_coverage = 0.9,
  max_target_seqs = 100,
  out_dir = tempfile("barcurateR_assign_"),
  blastn = "blastn",
  makeblastdb = "makeblastdb",
  table_name = "reference_final"
)

rb_build_species_matrix(
  assignment,
  asv_table,
  include_ambiguous = FALSE,
  unassigned = FALSE
)
```

## Arguments

- path:

  Path to a BLAST tabular output file or FASTA file.

- hits:

  A data frame of assignment hits (output from `rb_assign_blastn` or
  `rb_exact_hits`).

- min_identity:

  Minimum percent identity for a confident assignment.

- min_coverage:

  Minimum coverage threshold. In `rb_assign_blast()`, this is applied as
  a pre-filter on query coverage (`qcov`). In `rb_score_assignments()`,
  it is applied to combined coverage (`ccovs`) for high/low confidence
  assignment. Applied as a proportion (e.g., `0.9` for 90%). Combined
  coverage is calculated as the geometric mean of query and subject
  coverage.

- tie_bitscore:

  Tolerance for bitscore ties.

- asvs:

  A data frame of ASVs (output from `rb_read_fasta`).

- refs:

  A data frame of reference sequences.

- asv_fasta:

  Path to a FASTA file containing ASV sequences.

- out_dir:

  Directory for temporary BLASTn files.

- max_target_seqs:

  Maximum number of target sequences to retain per ASV in BLASTn.
  Defaults to 100 to ensure that ties across many congeneric species are
  not prematurely dropped by the BLAST algorithm.

- blastn:

  Path to the `blastn` executable.

- makeblastdb:

  Path to the `makeblastdb` executable.

- con:

  A DBI connection to a barcurateR SQLite database.

- db_path:

  Path to the SQLite database (used if `con` is NULL).

- marker:

  The marker gene to extract from the database (e.g., `"12S"`, `"COI"`).

- qc_flag:

  Character vector of QC flags to include from the reference database
  (e.g., `"pass"`, `"manual_review_needed"`). Passed to
  [`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md).
  Defaults to `"pass"`.

- method:

  Assignment method: `"blastn"` (uses BLAST+) or `"exact"`
  (dependency-free exact matching mode).

- table_name:

  Name of the reference table in the database (default
  `"reference_final"`).

- assignment:

  A data frame of scored assignments.

- asv_table:

  A data frame of ASV counts per sample.

- include_ambiguous:

  Logical. Include ambiguous and `"check_full_tie"` assignments in the
  matrix?

- unassigned:

  Logical. Include unassigned ASVs (`no_match`) in the matrix?

## Value

- `rb_assign_edna()`, `rb_assign_blastn()`, `rb_exact_hits()`, and
  `rb_score_assignments()` return a data frame of scored assignments.

- `rb_parse_blast_tabular()` and `rb_read_fasta()` return a parsed data
  frame.

- `rb_build_species_matrix()` returns a species-by-sample abundance
  matrix.

## Details

When `method = "blastn"`, BLAST+ must be installed and `blastn` and
`makeblastdb` must be available on the system path.

The assignment pipeline utilizes advanced alignment metrics to prevent
common BLAST biases. It calculates subject coverage (`scovs`), combined
coverage (`ccovs` = sqrt(qcov \* scovs)), and bitscore per base
(`bitscore_pb`). Hits are ranked using a robust hierarchy: `-bitscore`,
`-bitscore_pb`, `evalue`, `-ccovs`, `-pident`.

The `assignment_status` column classifies the result:

- `"unique"`: A single best species match.

- `"ambiguous"`: Multiple species tied, but with slight metric
  variations.

- `"check_full_tie"`: Multiple species tied with perfectly identical top
  metrics across the board, requiring manual biogeographical or
  taxonomic curation.

- `"no_match"`: No hits passed the coverage/identity filters.

For ties, the `tie_taxa` column contains a pipe-separated (`|`) string
of the tied species names (e.g., `"Species A|Species B"`). If there is
no tie, this column is `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assign ASVs using exact matching (no BLAST+ required)
assignment <- rb_assign_edna(
  asv_fasta = "asvs.fasta",
  db_path = "curated_ref.db",
  marker = "12S",
  method = "exact"
)

# Build a species-by-sample matrix
asv_counts <- read.csv("asv_table.csv")
species_matrix <- rb_build_species_matrix(assignment, asv_counts)
} # }
```
