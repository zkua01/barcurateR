# Sequence integrity checks for quality control

Functions to check individual sequences or groups of sequences for
common quality issues.

- `rb_check_codons()`: Checks coding sequences (COI) for stop codons and
  frameshifts.

- `rb_check_rrna_integrity()`: Checks rRNA sequences (12S, 16S) for
  abnormal length and the presence of any non-ACGT characters (including
  alignment gaps and ambiguous IUPAC codes like 'N').

- `rb_check_divergence()`: Detects phylogenetically divergent sequences
  within species/gene groups using MAFFT alignment and FastTree.

## Usage

``` r
rb_check_codons(sequence, gene, coding_genes = c("COI"), genetic_code = "2")

rb_check_rrna_integrity(
  sequence,
  gene,
  rrna_genes = c("12S", "16S"),
  min_length = 100
)

rb_check_divergence(
  data,
  species_col = "species",
  gene_col = "gene",
  sequence_col = "sequence",
  id_col = "unique_code",
  min_seqs = 5,
  outlier_mult = 5,
  outlier_floor = 0.02,
  mafft = "mafft",
  fasttree = "FastTree"
)
```

## Arguments

- sequence:

  A single DNA sequence string.

- gene:

  The gene/marker type (e.g., `"COI"`, `"12S"`, `"genome"`).

- coding_genes:

  Gene types to check for codon issues. Default `c("COI")`.

- genetic_code:

  Genetic code table (default `"2"` for vertebrate mtDNA).

- rrna_genes:

  Gene types to check for rRNA integrity.

- min_length:

  Minimum expected length for rRNA sequences.

- data:

  A data frame containing sequences and metadata.

- species_col, gene_col, sequence_col, id_col:

  Column names.

- min_seqs:

  Minimum sequences per group for divergence analysis.

- outlier_mult:

  Multiplier for median edge length to find outliers.

- outlier_floor:

  Minimum absolute edge length for outliers.

- mafft, fasttree:

  Paths to external executables.

## Value

- `rb_check_codons()` returns a list with `has_stop`, `frameshifted`,
  `best_frame`.

- `rb_check_rrna_integrity()` returns a list with `has_short`,
  `has_gaps`.

- `rb_check_divergence()` returns `data` with a `div_result` column.

## Details

Note: This function relies on shell pipes (\|) to stream MAFFT output
into FastTree. While this works on Unix-like systems, Windows users may
need to ensure their shell environment (e.g., Git Bash or WSL) is
properly configured, or install MSYS2, as the default cmd.exe may not
support these pipes.

## See also

Other quality control:
[`rb_qc_pipeline`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md),
[`rb_qc_screening`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md)
