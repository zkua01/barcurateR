# ============================================================
# R/phylogeny.R
# Reusable gene-tree building (alignment + distance + NJ tree +
# optional outgroup rooting). Deliberately does NOT include the
# Megalobrama-specific mislabel-detection logic — see project notes /
# examples/megalobrama_mislabel_correction.R for that, built on top of
# this function.
# ============================================================

# ------------------------------------------------------------
# rb_build_gene_tree()
#' Build a neighbor-joining tree for one gene/marker across the dataset
#'
#' @param outgroup Optional species name (or pattern, if
#'   `outgroup_exact = FALSE`) to root the tree on.
#' @param outgroup_exact If TRUE (default), `outgroup` is matched as an
#'   exact species via an anchored prefix match against tip labels
#'   ("^<outgroup>-"), avoiding accidental partial matches (e.g.
#'   outgroup = "sp2" no longer matches a tip for species "sp20"). Set
#'   FALSE for the old substring-anywhere behavior.
# ------------------------------------------------------------

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
        if (isTRUE(outgroup_exact)) "exact" else "substring", " match) — tree left unrooted.",
        call. = FALSE
      )
    }
  }

  list(tree = tree, method = "njs", n_tips = length(tree$tip.label))
}
