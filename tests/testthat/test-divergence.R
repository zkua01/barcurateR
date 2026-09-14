# ============================================================
# tests/testthat/test-divergence.R
#
# Tests for phylogenetic divergence checking.
#
# This module generalizes check_divergence_safe() from
# edna_ref_qc.R / edna_ref_qc_p2.R.
#
# Some tests require external tools:
#   - mafft
#   - FastTree
#
# Those tests will skip if the tools are not available.
# ============================================================


# ------------------------------------------------------------
# Helper: find external tool
# ------------------------------------------------------------

find_external_tool <- function(candidates) {
  for (tool in candidates) {
    path <- Sys.which(tool)
    if (nzchar(path)) {
      return(tool)
    }
  }
  ""
}


# ------------------------------------------------------------
# Helper: small test dataset
#
# Five similar sequences and one clearly divergent sequence.
# This is not a large real dataset, but it should be enough for
# mafft/FastTree to build a small tree locally.
# ------------------------------------------------------------

make_divergence_test_data <- function() {
  set.seed(42)
  
  seq_length <- 180
  
  base_seq <- paste(
    sample(c("A", "C", "G", "T"), seq_length, replace = TRUE),
    collapse = ""
  )
  
  mutate_sequence <- function(x, n_mutations = 3) {
    chars <- strsplit(x, "")[[1]]
    idx <- sample(seq_along(chars), n_mutations)
    chars[idx] <- sample(c("A", "C", "G", "T"), length(idx), replace = TRUE)
    paste(chars, collapse = "")
  }
  
  similar_sequences <- vapply(
    1:5,
    function(i) mutate_sequence(base_seq, n_mutations = 2),
    character(1)
  )
  
  divergent_sequence <- paste(
    sample(c("A", "C", "G", "T"), seq_length, replace = TRUE),
    collapse = ""
  )
  
  data.frame(
    species = "SpeciesA",
    gene = "COI",
    unique_code = c(paste0("similar_", 1:5), "divergent_1"),
    sequence = c(similar_sequences, divergent_sequence),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Offline tests: these should not require mafft/FastTree
# ------------------------------------------------------------

test_that("rb_check_divergence handles groups with insufficient sequences", {
  
  small_data <- make_divergence_test_data()[1:2, ]
  
  result <- rb_check_divergence(
    small_data,
    species_col = "species",
    gene_col = "gene",
    sequence_col = "sequence",
    id_col = "unique_code",
    min_seqs = 5
  )
  
  expect_s3_class(result, "data.frame")
  
  # Depending on implementation, the function may return zero rows
  # for groups below the threshold, or it may return rows labelled
  # insufficient_data / not_checked.
  if (nrow(result) > 0) {
    expect_true(
      all(result$div_result %in% c("insufficient_data", "not_checked"))
    )
  }
})


test_that("rb_check_divergence detects extreme length variation", {
  
  length_data <- data.frame(
    species = "SpeciesA",
    gene = "COI",
    unique_code = paste0("seq_", 1:5),
    sequence = c(
      strrep("A", 20),
      strrep("A", 20),
      strrep("A", 20),
      strrep("A", 20),
      strrep("A", 4000)
    ),
    stringsAsFactors = FALSE
  )
  
  result <- rb_check_divergence(
    length_data,
    species_col = "species",
    gene_col = "gene",
    sequence_col = "sequence",
    id_col = "unique_code",
    min_seqs = 2
  )
  
  expect_s3_class(result, "data.frame")
  
  expect_true("div_result" %in% names(result))
  
  expect_true(
    all(result$div_result == "high_length_variation")
  )
})


# ------------------------------------------------------------
# External tool tests: require mafft and FastTree
# ------------------------------------------------------------

test_that("rb_check_divergence runs mafft/FastTree and returns valid statuses", {
  
  mafft_bin <- find_external_tool(c("mafft"))
  fasttree_bin <- find_external_tool(c("FastTree", "fasttree"))
  
  skip_if_not(nzchar(mafft_bin), "mafft is not installed")
  skip_if_not(nzchar(fasttree_bin), "FastTree is not installed")
  skip_if_not_installed("ape")
  skip_on_cran()
  
  test_data <- make_divergence_test_data()
  
  result <- rb_check_divergence(
    test_data,
    species_col = "species",
    gene_col = "gene",
    sequence_col = "sequence",
    id_col = "unique_code",
    min_seqs = 5,
    mafft = mafft_bin,
    fasttree = fasttree_bin
  )
  
  expect_s3_class(result, "data.frame")
  
  expect_true(all(c("unique_code", "div_result") %in% names(result)))
  
  # All input sequences should be represented in the output.
  expect_true(all(test_data$unique_code %in% result$unique_code))
  
  valid_statuses <- c(
    "pass",
    "divergent",
    "tree_error",
    "processing_error",
    "high_length_variation",
    "insufficient_data",
    "not_checked"
  )
  
  expect_true(all(result$div_result %in% valid_statuses))
  
  # On a healthy local installation, this small test should not produce
  # processing errors or tree errors.
  expect_false(any(result$div_result == "processing_error"))
  expect_false(any(result$div_result == "tree_error"))
})


test_that("rb_check_divergence can flag an obvious divergent sequence", {
  
  mafft_bin <- find_external_tool(c("mafft"))
  fasttree_bin <- find_external_tool(c("FastTree", "fasttree"))
  
  skip_if_not(nzchar(mafft_bin), "mafft is not installed")
  skip_if_not(nzchar(fasttree_bin), "FastTree is not installed")
  skip_if_not_installed("ape")
  skip_on_cran()
  
  test_data <- make_divergence_test_data()
  
  result <- rb_check_divergence(
    test_data,
    species_col = "species",
    gene_col = "gene",
    sequence_col = "sequence",
    id_col = "unique_code",
    min_seqs = 5,
    mafft = mafft_bin,
    fasttree = fasttree_bin
  )
  
  # If tree building succeeded, we expect the obviously divergent
  # sequence to be flagged.
  if (!any(result$div_result %in% c("tree_error", "processing_error"))) {
    
    expect_true("divergent" %in% result$div_result)
    
    divergent_status <- result$div_result[result$unique_code == "divergent_1"]
    
    expect_equal(divergent_status, "divergent")
  }
})


test_that("rb_check_divergence handles multiple species/gene groups", {
  
  mafft_bin <- find_external_tool(c("mafft"))
  fasttree_bin <- find_external_tool(c("FastTree", "fasttree"))
  
  skip_if_not(nzchar(mafft_bin), "mafft is not installed")
  skip_if_not(nzchar(fasttree_bin), "FastTree is not installed")
  skip_if_not_installed("ape")
  skip_on_cran()
  
  group_a <- make_divergence_test_data()
  
  group_b <- make_divergence_test_data()
  group_b$species <- "SpeciesB"
  group_b$unique_code <- paste0("B_", group_b$unique_code)
  
  # Make SpeciesB too small to be checked.
  group_b <- group_b[1:2, ]
  
  combined <- rbind(group_a, group_b)
  
  result <- rb_check_divergence(
    combined,
    species_col = "species",
    gene_col = "gene",
    sequence_col = "sequence",
    id_col = "unique_code",
    min_seqs = 5,
    mafft = mafft_bin,
    fasttree = fasttree_bin
  )
  
  expect_s3_class(result, "data.frame")
  
  # SpeciesA should be checked.
  expect_true(any(result$unique_code %in% group_a$unique_code))
  
  # SpeciesB may either be absent from the result or labelled as
  # insufficient_data / not_checked.
  if (any(result$unique_code %in% group_b$unique_code)) {
    b_rows <- result[result$unique_code %in% group_b$unique_code, ]
    expect_true(
      all(b_rows$div_result %in% c("insufficient_data", "not_checked"))
    )
  }
})


# ------------------------------------------------------------
# Graceful failure test: invalid external tool paths
# ------------------------------------------------------------

test_that("rb_check_divergence handles missing external tools gracefully", {
  
  skip_if_not_installed("ape")
  
  test_data <- make_divergence_test_data()
  
  result <- tryCatch(
    rb_check_divergence(
      test_data,
      species_col = "species",
      gene_col = "gene",
      sequence_col = "sequence",
      id_col = "unique_code",
      min_seqs = 5,
      mafft = "not_a_real_mafft_binary",
      fasttree = "not_a_real_fasttree_binary"
    ),
    error = function(e) e
  )
  
  # Either the function errors informatively, or it returns
  # tree_error / processing_error rows.
  if (inherits(result, "error")) {
    expect_match(
      conditionMessage(result),
      "mafft|FastTree|fasttree|not found|executable",
      ignore.case = TRUE
    )
  } else {
    expect_s3_class(result, "data.frame")
    expect_true(
      all(result$div_result %in% c("tree_error", "processing_error"))
    )
  }
})