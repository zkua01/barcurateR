# Prepare sequences for de novo alignment (internal)

Removes alignment padding and whitespace, and converts IUPAC ambiguity
codes to `N`, for algorithms that require an ungapped ACGTN input. Never
used to modify stored sequences.

## Usage

``` r
rb_prepare_for_alignment(sequence)
```

## Arguments

- sequence:

  Character vector.

## Value

Character vector of the same length as `sequence`, uppercased, with
alignment padding and whitespace removed and all non-ACGT characters
replaced by `N`.

## Details

Intended solely for the input boundary of de novo alignment, where an
ungapped ACGTN alphabet is required. **Do not apply to already-aligned
sequences**: gap characters carry the column correspondence between
taxa, and removing them invalidates per-column operations such as
[`Biostrings::consensusMatrix()`](https://rdrr.io/pkg/Biostrings/man/letterFrequency.html)
and equal-width distance calculations.
