# Prepare and trim raw DNA sequences

Utility functions for preparing raw DNA sequences during source
standardization.

- `rb_reverse_complement()`: Returns the reverse complement of a DNA
  sequence.

- `rb_remove_primers()`: Removes matching forward and reverse primer
  sequences from the ends of a sequence, if present.

## Usage

``` r
rb_reverse_complement(dna_seq)

rb_remove_primers(
  sequence,
  f_primer_seqs = NA_character_,
  r_primer_seqs = NA_character_
)
```

## Arguments

- dna_seq:

  Character vector of DNA sequences.

- sequence:

  Character vector of DNA sequences to trim.

- f_primer_seqs:

  Character vector of forward primer sequences. Use `NA` if no forward
  primers should be removed.

- r_primer_seqs:

  Character vector of reverse primer sequences. Use `NA` if no reverse
  primers should be removed.

## Value

A character vector of transformed sequences.

## Details

`rb_remove_primers()` performs simple exact matching at sequence ends.
Forward primers are removed from the start of the sequence. Reverse
primers are reverse-complemented and removed from the end of the
sequence.

## See also

Other parsing:
[`rb_ambiguity_utils`](https://zkua01.github.io/barcurateR/reference/rb_ambiguity_utils.md),
[`rb_combine_resolve`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md),
[`rb_parse_sources`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rb_reverse_complement("ATGC")

rb_remove_primers(
  "ACGTACGTACGT",
  f_primer_seqs = "ACGT",
  r_primer_seqs = NA
)
} # }
```
