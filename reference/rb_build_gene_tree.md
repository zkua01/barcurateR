# Build a neighbor-joining gene tree for a marker

Builds a neighbor-joining phylogenetic tree for a single gene or marker
across all sequences in the dataset. Sequences are aligned with
[`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html),
pairwise distances are computed with
[`ape::dist.dna()`](https://rdrr.io/pkg/ape/man/dist.dna.html), and the
tree is constructed with
[`ape::njs()`](https://rdrr.io/pkg/ape/man/njs.html).

Optionally, the tree can be rooted on a user-specified outgroup.

## Usage

``` r
rb_build_gene_tree(
  data,
  gene_type,
  gene_col = "seq_type",
  species_col = "species",
  sequence_col = "sequence",
  id_col = "unique_code",
  outgroup = NULL,
  outgroup_exact = TRUE,
  dna_model = "K80",
  min_seqs = 3,
  min_length = 100
)
```

## Arguments

- data:

  A standardized data frame containing at least the columns specified by
  `gene_col`, `species_col`, `sequence_col`, and `id_col`.

- gene_type:

  The marker/gene to build the tree for (e.g., `"COI"`, `"12S"`).

- gene_col:

  Name of the column containing marker/gene labels. Default
  `"seq_type"`.

- species_col:

  Name of the species column. Default `"species"`.

- sequence_col:

  Name of the sequence column. Default `"sequence"`.

- id_col:

  Name of the unique sequence identifier column. Default
  `"unique_code"`.

- outgroup:

  Optional species name to root the tree on. If `NULL`, the tree is left
  unrooted.

- outgroup_exact:

  Logical. If `TRUE` (default), `outgroup` is matched as an anchored
  prefix against tip labels (`"^<outgroup>-"`), avoiding accidental
  partial matches. Set `FALSE` for substring-anywhere matching.

- dna_model:

  DNA substitution model passed to
  [`ape::dist.dna()`](https://rdrr.io/pkg/ape/man/dist.dna.html).
  Default `"K80"` (Kimura 2-parameter).

- min_seqs:

  Minimum number of sequences required to build a tree. Default `3`.

- min_length:

  Minimum sequence length (in bp) after gap removal. Default `100`.

## Value

A list containing:

- tree:

  A rooted or unrooted `phylo` object.

- method:

  The tree-building method used (`"njs"`).

- n_tips:

  The number of tips in the tree.

Returns `NULL` (with a warning) if fewer than `min_seqs` sequences
remain after filtering. Alignment or distance computation errors will
propagate to the caller.

## Details

The function performs the following steps:

1.  Filters `data` to the requested `gene_type`.

2.  Removes duplicate sequences.

3.  Strips alignment gap characters (`"-"`).

4.  Filters sequences shorter than `min_length`.

5.  Aligns sequences using
    [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html).

6.  Computes pairwise distances with
    [`ape::dist.dna()`](https://rdrr.io/pkg/ape/man/dist.dna.html).

7.  Builds a neighbor-joining tree with
    [`ape::njs()`](https://rdrr.io/pkg/ape/man/njs.html).

8.  Optionally roots the tree on the specified outgroup.

Tip labels are formatted as `"<species>-<unique_code>"`, with spaces in
species names replaced by hyphens.

This function deliberately does not include taxon-specific mislabel
detection logic. For that, build on top of this function using
project-specific heuristics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build an unrooted COI tree
tree_result <- rb_build_gene_tree(final_data, gene_type = "COI")
plot(tree_result$tree)

# Build a rooted 12S tree using an outgroup
tree_result <- rb_build_gene_tree(
  final_data,
  gene_type = "12S",
  outgroup = "Alligator sinensis"
)
plot(tree_result$tree)
} # }
```
