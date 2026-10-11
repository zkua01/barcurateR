# Export curated reference sequences to common metabarcoding formats

Functions to export curated reference sequences and taxonomy to
FASTA-based formats used by common metabarcoding pipelines.

- `rb_export_fasta()`: General dispatcher that writes to one of the
  supported formats.

- `rb_export_blastn()`: Writes BLAST-style FASTA headers.

- `rb_export_dada2()`: Writes DADA2-style FASTA headers.

- `rb_export_qiime2()`: Writes a QIIME2-compatible FASTA file and a
  separate taxonomy TSV file.

- `rb_export_usearch()`: Writes USEARCH/VSEARCH-style FASTA headers.

## Usage

``` r
rb_export_fasta(
  data,
  file,
  format = c("blastn", "dada2", "qiime2", "usearch"),
  taxonomy_file = NULL
)

rb_export_blastn(data, fasta_file)

rb_export_dada2(data, fasta_file)

rb_export_qiime2(data, fasta_file, taxonomy_file)

rb_export_usearch(data, fasta_file)
```

## Arguments

- data:

  A data frame containing at least `sequence` and the seven taxonomic
  ranks: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, and
  `species`.

- file:

  Output file path. Used by `rb_export_fasta()`.

- format:

  Export format. One of `"blastn"`, `"dada2"`, `"qiime2"`, or
  `"usearch"`. Used by `rb_export_fasta()`.

- taxonomy_file:

  Output taxonomy TSV file path. Required for `rb_export_qiime2()` and
  `rb_export_fasta(format = "qiime2")`.

- fasta_file:

  Output FASTA file path. Used by the format-specific export functions.

## Value

- `rb_export_fasta()`, `rb_export_blastn()`, `rb_export_dada2()`, and
  `rb_export_usearch()` return the output FASTA path invisibly.

- `rb_export_qiime2()` returns an invisible list containing `fasta` and
  `taxonomy` file paths.

## Details

All export functions require full taxonomy from kingdom to species. They
do not perform QC filtering themselves.

To export only QC-passing sequences, filter the data first using
[`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md)
or
[`rb_filter_reference()`](https://zkua01.github.io/barcurateR/reference/rb_query.md).
For example:

    pass_refs <- rb_get_sequences(con, qc_flag = "pass")
    rb_export_dada2(pass_refs, "reference_dada2.fasta")

Header styles differ by format:

|         |                                                 |
|---------|-------------------------------------------------|
| Format  | Example header style                            |
| BLAST   | `>seq_id;Kingdom;Phylum;Class;...`              |
| DADA2   | `>seq_id;tax=Kingdom;Phylum;Class;...`          |
| USEARCH | `>seq_id tax=k__Kingdom;p__Phylum;c__Class;...` |
| QIIME2  | Plain FASTA headers plus separate taxonomy TSV  |

## Examples

``` r
if (FALSE) { # \dontrun{
con <- rb_connect()
refs <- rb_get_sequences(con, qc_flag = "pass")
rb_disconnect(con)

rb_export_dada2(refs, "ref_dada2.fasta")
rb_export_qiime2(refs, "ref_qiime2.fasta", "ref_taxonomy.tsv")
} # }
```
