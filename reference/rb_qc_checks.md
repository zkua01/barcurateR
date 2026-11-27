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
rb_check_codons(
  sequence,
  gene,
  coding_genes = c("COI", "CYTB", "complete_genome", "genome", "multi_marker"),
  genetic_code = "2",
  max_internal_stops = 2L
)

rb_check_rrna_integrity(
  sequence,
  gene,
  rrna_genes = c("12S", "16S"),
  min_length = 100L,
  max_ambiguous_frac = NULL
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

- max_internal_stops:

  Integer scalar (default `2L`). Maximum number of internal stop codons
  tolerated in the best of the six reading frames before the record is
  flagged `frameshifted`. The terminal stop of a complete CDS is not
  counted. `has_stop` is set independently whenever any internal stop is
  present, so lowering this value only tightens the `frameshifted` flag,
  not `internal_stop`.

- rrna_genes:

  Gene types to check for rRNA integrity.

- min_length:

  Minimum expected length for rRNA sequences.

- max_ambiguous_frac:

  Numeric scalar in `[0, 1]`, or `NULL` (default). Optional exclusion
  threshold on the proportion of ambiguous (IUPAC) bases among called
  bases (`frac_ambig`). When supplied, records with
  `frac_ambig > max_ambiguous_frac` receive the `ambiguous_excess` flag
  and are excluded from the passing set. When `NULL`, ambiguity never
  excludes a record; it is reported only through the informational
  `ambiguous_bases` note and the `n_ambiguous` / `frac_ambig` metadata
  columns.

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

- `rb_check_codons()` returns a `data.frame` with one row per input
  sequence:

  - has_stop:

    Logical. `TRUE` if the best reading frame contains at least one
    internal stop codon. `NA` for non-coding markers.

  - frameshifted:

    Logical. `TRUE` iif `n_internal_stops` exceeds `max_internal_stops`.
    `NA` for non-coding markers.

  - n_internal_stops:

    Integer. Internal stop count in the best frame, excluding a terminal
    stop. `NA` for non-coding markers.

  - best_frame:

    Integer. Index of the selected frame among those evaluated (see
    Details). `NA` for non-coding markers.

  - coding_strand:

    Character. `"F"` if the best frame is on the supplied strand, `"R"`
    if on its reverse complement. `NA` for non-coding markers.

  - frame_complete:

    Logical. Whether the cleaned sequence length is a multiple of three
    (informational only). `NA` for non-coding markers.

  - codon_check_run:

    Logical tri-state

- `rb_check_rrna_integrity()` returns a `data.frame` with one row per
  input sequence:

  - has_gaps:

    Logical. `TRUE` if any gap character or IUPAC code is present.
    Informational; does not exclude a record. `NA` for non-rRNA markers.

  - has_short:

    Logical. `TRUE` if `n_base` is below `min_length`. `NA` for non-rRNA
    markers.

  - n_gap_chars:

    Integer. Count of `-` and `.` characters.

  - n_ambiguous:

    Integer. Count of IUPAC ambiguity codes.

  - n_base:

    Integer. Biological length, `n_acgt + n_ambiguous`.

  - frac_ambig:

    Numeric. `n_ambiguous / n_base`; `NA` when `n_base` is 0.

  - max_gap_run:

    Integer. Length of the longest contiguous gap run.

- `rb_check_divergence()` returns `data` with a `div_result` column.

## Details

Note: This function relies on shell pipes (\|) to stream MAFFT output
into FastTree. While this works on Unix-like systems, Windows users may
need to ensure their shell environment (e.g., Git Bash or WSL) is
properly configured, or install MSYS2, as the default cmd.exe may not
support these pipes.

Translation is attempted in six frames: three forward offsets and three
on the reverse complement. The frame with the fewest internal stops is
selected, so records deposited in reverse orientation are handled
correctly. Frames shorter than three bases are skipped; for very short
inputs `best_frame` therefore indexes the reduced set, and
`coding_strand` is the more reliable identifier.

Alignment padding (`-`, `.`) and whitespace are removed before
translation; IUPAC ambiguity codes are retained and resolved with
`Biostrings::translate(if.fuzzy.codon = "solve")`, so an ambiguous codon
does not abort the check. A stop in the final codon is treated as a
legitimate termination and excluded from the internal count.

Stops are those of `genetic_code`, which defaults to `"2"` (vertebrate
mitochondrial).

Gap characters and ambiguity codes are reported separately because they
mean different things. `-` and `.` may be alignment padding, a genuine
deletion relative to other sequences, or undetermined bases; these
cannot be distinguished from the character alone, so nothing is stripped
or altered here. IUPAC codes represent real base-calling uncertainty and
are retained.

Neither condition excludes a record by default. Presence of gap or
ambiguous characters is recorded in `has_gaps` and surfaced as the
informational `sequence_gap` / `ambiguous_bases` notes by
[`rb_compile_qc_flags()`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md).
Only `has_short` is excluding, and it is assessed on biological length
so that alignment padding cannot inflate a short record past the
threshold.

## See also

Other quality control:
[`rb_qc_pipeline`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md),
[`rb_qc_screening`](https://zkua01.github.io/barcurateR/reference/rb_qc_screening.md)

## Examples

``` r
# Non-coding markers: check not applicable
rb_check_codons("ACGTACGTACGT", gene = "12S")
#>   has_stop frameshifted n_internal_stops best_frame coding_strand
#> 1       NA           NA               NA         NA          <NA>
#>   frame_complete codon_check_run
#> 1             NA              NA

# Coding marker: check completes, length is a multiple of three
res <- rb_check_codons(strrep("ACGT", 30), gene = "COI")
res$codon_check_run
#> [1] TRUE
res$frame_complete
#> [1] TRUE
res$coding_strand
#> [1] "F"

#' # Clean 200 bp reference
rb_check_rrna_integrity(strrep("ACGT", 50), gene = "12S")
#>   has_gaps has_short n_gap_chars n_ambiguous n_base frac_ambig max_gap_run
#> 1    FALSE     FALSE           0           0    200          0           0

# One ambiguity code: informational, still passing
s <- strrep("ACGT", 50)
s <- paste0(substr(s, 1, 100), "S", substr(s, 101, 200))
rb_check_rrna_integrity(s, gene = "12S")
#>   has_gaps has_short n_gap_chars n_ambiguous n_base  frac_ambig max_gap_run
#> 1     TRUE     FALSE           0           1    201 0.004975124           0

# Padding does not count toward biological length
padded <- paste0(strrep("A", 90), strrep("-", 50))
rb_check_rrna_integrity(padded, gene = "12S")$has_short   # TRUE
#> [1] TRUE
```
