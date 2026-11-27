# Characterise nucleotide composition of DNA sequences

Reports the character-class composition of sequences without modifying
them. Separates unambiguous bases, IUPAC ambiguity codes, gap characters
and any other characters, so that quality control can distinguish
formatting artefacts from genuine base-calling uncertainty.

## Usage

``` r
rb_sequence_composition(sequence)
```

## Arguments

- sequence:

  Character vector of nucleotide sequences. Case is ignored. `NA` values
  are treated as empty strings.

## Value

A [data.frame](https://rdrr.io/r/base/data.frame.html) with one row per
input sequence:

- n_char:

  Total characters in the stored string.

- n_acgt:

  Count of unambiguous bases (A, C, G, T).

- n_dash:

  Count of hyphen (`-`) characters.

- n_dot:

  Count of period (`.`) characters.

- n_ambig:

  Count of IUPAC ambiguity codes (R, Y, K, M, S, W, B, D, H, V, N).

- n_other:

  Count of characters belonging to none of the above classes.

- n_gap_runs:

  Number of contiguous runs of gap characters.

- max_gap_run:

  Length of the longest contiguous gap run; 0 if none.

- n_base:

  Biological length, `n_acgt + n_ambig`. Excludes gap characters.

- frac_ambig:

  `n_ambig / n_base`; `NA` when `n_base` is 0.

## Details

Gap characters are **not** stripped or altered. Their meaning is
context-dependent: they may be alignment padding, a genuine deletion
relative to other sequences in an alignment, or undetermined bases in
the source record. These cases cannot be distinguished from the
character alone, so they are reported separately and any normalisation
is left to an explicit, user-initiated step.

IUPAC codes are retained as information. Use `n_base` as the denominator
for ambiguity fractions and as the length measure for marker-level
checks; `n_char` reflects the stored string and may be inflated by
padding.

## See also

Other parsing:
[`rb_ambiguity_utils`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md),
[`rb_combine_resolve`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md),
[`rb_parse_sources`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md),
[`rb_sequence_prep`](https://zkua01.github.io/barcurateR/reference/rb_sequence_prep.md)

## Examples

``` r
x <- rb_sequence_composition(c("ACGTACGT", "ACGT-ACGT", "ACGTNACGT", "ACGT..ACGT"))
x$n_char   # 8  9  9 10  (stored string)
#> [1]  8  9  9 10
x$n_base   # 8  8  9  8  (biological length)
#> [1] 8 8 9 8
```
