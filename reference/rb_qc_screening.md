# Screen reference sequences for contaminants and NUMTs

Functions to identify sequences that match known contaminants or nuclear
mitochondrial pseudogenes (NUMTs).

- `rb_screen_contaminants()`: Screens sequences against a BLAST
  database.

- `rb_screen_numts()`: Screens sequences against a FASTA file of known
  NUMTs using exact substring matching.

## Usage

``` r
rb_screen_contaminants(
  sequences,
  blast_db,
  perc_identity = 95,
  evalue = 1e-20,
  min_pident = 95,
  min_length = 100,
  threads = 1,
  blastn = "blastn"
)

rb_screen_numts(sequences, numt_fasta, min_match_length = 50)
```

## Arguments

- sequences:

  Character vector of DNA sequences to screen.

- blast_db:

  Path prefix to a BLAST nucleotide database.

- perc_identity:

  Minimum percent identity passed directly to BLAST -perc_identify flag
  as a hard pre-filter for BLAST hits.

- evalue:

  E-value threshold for BLAST hits.

- min_pident:

  Post-filter minimum percent identity threshold applied to the BLAST
  results to flag as contaminant.

- min_length:

  Minimum alignment length to flag as contaminant.

- threads:

  Number of BLAST threads.

- blastn:

  Path to the `blastn` executable.

- numt_fasta:

  Path to a FASTA file of known NUMT sequences.

- min_match_length:

  Sequences shorter than this are not screened.

## Value

A logical vector where `TRUE` indicates the sequence was flagged.

## See also

Other quality control:
[`rb_qc_checks`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md),
[`rb_qc_pipeline`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md)
