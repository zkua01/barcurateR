# barcurateR: Standardize, Curate, Audit, and Deploy DNA Reference Databases

Provides tools to standardize, curate, query, audit, and export DNA
reference databases, and assign amplicon sequence variants for
ecological metabarcoding in a reproducible manner.

## Main workflow

1.  Parse and standardize source sequences
    ([`rb_parse_source_table()`](https://zkua01.github.io/barcurateR/reference/rb_parse_sources.md))

2.  Combine sources and resolve ambiguities
    ([`rb_combine_sources()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md),
    [`rb_resolve_ambiguous()`](https://zkua01.github.io/barcurateR/reference/rb_combine_resolve.md))

3.  Run the full curation pipeline
    ([`rb_curate_reference()`](https://zkua01.github.io/barcurateR/reference/rb_curate_reference.md))

4.  Query and export the curated database
    ([`rb_get_sequences()`](https://zkua01.github.io/barcurateR/reference/rb_query.md),
    [`rb_export_fasta()`](https://zkua01.github.io/barcurateR/reference/rb_export.md))

5.  Assign eDNA sequences
    ([`rb_assign_edna()`](https://zkua01.github.io/barcurateR/reference/rb_assignment.md))

## Vignette

Please see the vignettes for a detailed walkthrough of the curation and
assignment workflows (`browseVignettes("barcurateR")`)

## See also

Useful links:

- <https://zkua01.github.io/barcurateR/>

- <https://github.com/zkua01/barcurateR>

- Report bugs at <https://github.com/zkua01/barcurateR/issues>

## Author

**Maintainer**: Zi Xun Kua <zixunkua@gmail.com>

Authors:

- Zi Xun Kua <zixunkua@gmail.com>

- Katja Pogačar

- Guohuan Su
