# ============================================================
#' Build a neighbor-joining gene tree for a marker
#'
#' @description
#' Builds a neighbor-joining phylogenetic tree for a single gene or marker
#' across all sequences in the dataset. Sequences are aligned with
#' `DECIPHER::AlignSeqs()`, pairwise distances are computed with
#' `ape::dist.dna()`, and the tree is constructed with `ape::njs()`.
#'
#' Optionally, the tree can be rooted on a user-specified outgroup.
#'
#' @param data A standardized data frame containing at least the columns
#'   specified by `gene_col`, `species_col`, `sequence_col`, and `id_col`.
#' @param gene_type The marker/gene to build the tree for (e.g., `"COI"`,
#'   `"12S"`).
#' @param gene_col Name of the column containing marker/gene labels.
#'   Default `"seq_type"`.
#' @param species_col Name of the species column. Default `"species"`.
#' @param sequence_col Name of the sequence column. Default `"sequence"`.
#' @param id_col Name of the unique sequence identifier column.
#'   Default `"unique_code"`.
#' @param outgroup Optional species name to root the tree on. If `NULL`,
#'   the tree is left unrooted.
#' @param outgroup_exact Logical. If `TRUE` (default), `outgroup` is matched
#'   as an anchored prefix against tip labels (`"^<outgroup>-"`), avoiding
#'   accidental partial matches. Set `FALSE` for substring-anywhere matching.
#' @param dna_model DNA substitution model passed to `ape::dist.dna()`.
#'   Default `"K80"` (Kimura 2-parameter).
#' @param min_seqs Minimum number of sequences required to build a tree.
#'   Default `3`.
#' @param min_length Minimum sequence length (in bp) after gap removal.
#'   Default `100`.
#'
#' @return A list containing:
#' \describe{
#'   \item{tree}{A rooted or unrooted `phylo` object.}
#'   \item{method}{The tree-building method used (`"njs"`).}
#'   \item{n_tips}{The number of tips in the tree.}
#' }
#' Returns `NULL` (with a warning) if fewer than `min_seqs` sequences
#' remain after filtering.
#'
#' @details
#' The function performs the following steps:
#'
#' 1. Filters `data` to the requested `gene_type`.
#' 2. Removes duplicate sequences.
#' 3. Strips alignment gap characters (`"-"`).
#' 4. Filters sequences shorter than `min_length`.
#' 5. Aligns sequences using `DECIPHER::AlignSeqs()`.
#' 6. Computes pairwise distances with `ape::dist.dna()`.
#' 7. Builds a neighbor-joining tree with `ape::njs()`.
#' 8. Optionally roots the tree on the specified outgroup.
#'
#' Tip labels are formatted as `"<species>-<unique_code>"`, with spaces
#' in species names replaced by hyphens.
#'
#' This function deliberately does not include taxon-specific mislabel
#' detection logic. For that, build on top of this function using
#' project-specific heuristics.
#'
#' @examples
#' \dontrun{
#' # Build an unrooted COI tree
#' tree_result <- rb_build_gene_tree(final_data, gene_type = "COI")
#' plot(tree_result$tree)
#'
#' # Build a rooted 12S tree using an outgroup
#' tree_result <- rb_build_gene_tree(
#'   final_data,
#'   gene_type = "12S",
#'   outgroup = "Alligator sinensis"
#' )
#' plot(tree_result$tree)
#' }
#'
#' @family phylogeny
#' @export
rb_build_gene_tree <- function(data, gene_type, gene_col = "seq_type",
                                species_col = "species", sequence_col = "sequence",
                                id_col = "unique_code", outgroup = NULL,
                                outgroup_exact = TRUE,
                                dna_model = "K80", min_seqs = 3, min_length = 100) {
  rb_required_columns(data, c(gene_col, species_col, sequence_col, id_col))

  gene_data <- data[data[[gene_col]] == gene_type, , drop = FALSE]
  gene_data <- gene_data[!duplicated(gene_data[[sequence_col]]), , drop = FALSE]
  gene_data[[sequence_col]] <- gsub("-", "", gene_data[[sequence_col]])
  gene_data <- gene_data[nchar(gene_data[[sequence_col]]) >= min_length, , drop = FALSE]

  if (nrow(gene_data) < min_seqs) {
    warning("Not enough sequences for ", gene_type, call. = FALSE)
    return(NULL)
  }

  tree_code <- paste0(gsub(" ", "-", gene_data[[species_col]]), "-", gene_data[[id_col]])
  dna_seqs <- Biostrings::DNAStringSet(stats::setNames(gene_data[[sequence_col]], tree_code))

  alignment <- DECIPHER::AlignSeqs(dna_seqs, iteration = 3, refinement = 3, verbose = FALSE)
  alignment_bin <- ape::as.DNAbin(alignment)

  dist_mat <- ape::dist.dna(alignment_bin, model = dna_model, pairwise.deletion = FALSE, as.matrix = TRUE)
  if (any(is.na(dist_mat))) {
    max_dist <- max(dist_mat, na.rm = TRUE)
    dist_mat[is.na(dist_mat)] <- max_dist * 1.1
  }

  tree <- ape::njs(dist_mat)

  if (!is.null(outgroup)) {
    outgroup_tips <- if (isTRUE(outgroup_exact)) {
      prefix <- paste0("^", gsub(" ", "-", outgroup), "-")
      grep(prefix, tree$tip.label, value = TRUE)
    } else {
      grep(outgroup, tree$tip.label, value = TRUE)
    }
    if (length(outgroup_tips) > 0) {
      tree <- ape::root(tree, outgroup = outgroup_tips[1])
    } else {
      warning(
        "Outgroup pattern '", outgroup, "' not found in tree tips (",
        if (isTRUE(outgroup_exact)) "exact" else "substring", " match) -- tree left unrooted.",
        call. = FALSE
      )
    }
  }

  list(tree = tree, method = "njs", n_tips = length(tree$tip.label))
}
