# Changelog

## barcurateR 0.3.0

### Key changes

- `sequence_gap` is now informational: gap characters and IUPAC codes no
  longer exclude a record from the passing set. Re-curate existing
  databases.
- `frameshifted` is now determined by internal stop codons in the best
  of six reading frames, not by sequence length modulo 3. Partial
  in-frame CDS (including the 658 bp Folmer barcode) are no longer
  flagged.
- [`rb_compile_qc_flags()`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md)
  replaces the `extra_flag_fn` mechanism. `qc_flag` contains only
  excluding flags; informational conditions are in the new `qc_note`
  column.

### New features

- [`rb_sequence_composition()`](https://zkua01.github.io/barcurateR/reference/rb_sequence_composition.md)
  reports character-class composition without modifying sequences.
- [`rb_check_codons()`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md)
  evaluates six frames, tolerates ambiguous codons, and ignores terminal
  stops. New metadata: `n_internal_stops`, `best_frame`,
  `coding_strand`, `frame_complete`, `codon_check_run`.
- [`rb_check_rrna_integrity()`](https://zkua01.github.io/barcurateR/reference/rb_qc_checks.md)
  reports `n_gap_chars`, `n_ambiguous`, `n_base`, `frac_ambig`,
  `max_gap_run`; length checks use biological length.
- Optional `max_ambiguous_frac` exclusion threshold in
  [`rb_compile_qc_flags()`](https://zkua01.github.io/barcurateR/reference/rb_qc_pipeline.md).

### Bug fixes

- Alignment padding no longer crashes
  [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html)
  in
  [`rb_build_gene_tree()`](https://zkua01.github.io/barcurateR/reference/rb_build_gene_tree.md)
  and
  [`rb_run_barcode_gap()`](https://zkua01.github.io/barcurateR/reference/rb_barcode_gap_analysis.md).
- Classifier length feature uses biological rather than stored string
  length.
