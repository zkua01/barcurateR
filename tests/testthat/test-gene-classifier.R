# ============================================================
# tests/testthat/test-gene-classifier.R
#
# Tests for the ML gene classifier and the QC round-trip.
# ============================================================


# ------------------------------------------------------------
# rb_sequence_features(): feature extraction
# ------------------------------------------------------------

test_that("rb_sequence_features computes correct GC/AT content and skew", {
  
  # Simple known case: 50% GC, 50% AT
  result <- rb_sequence_features("AATTGGCC")
  
  expect_s3_class(result, "data.frame")
  expect_equal(result$length, 8)
  expect_equal(result$gc_content, 0.5)
  expect_equal(result$at_content, 0.5)
  expect_equal(result$gc_skew, 0)    # equal G and C
  expect_equal(result$at_skew, 0)    # equal A and T
})


test_that("rb_sequence_features handles vectorized input", {
  
  result <- rb_sequence_features(c("AATTGGCC", "GGGGCCCC"))
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(result$length, c(8, 8))
  expect_equal(result$gc_content, c(0.5, 1.0))
  expect_equal(result$at_content, c(0.5, 0.0))
})


test_that("rb_sequence_features handles empty sequence without error", {
  
  result <- rb_sequence_features("")
  
  expect_s3_class(result, "data.frame")
  expect_equal(result$length, 0)
  expect_equal(result$gc_content, 0)
  expect_equal(result$at_content, 0)
})


test_that("rb_sequence_features computes skew correctly", {
  
  # GGGGAAAA: G=4, C=0, A=4, T=0
  result <- rb_sequence_features("GGGGAAAA")
  
  expect_equal(result$gc_skew, 1)    # (G-C)/(G+C) = (4-0)/(4+0) = 1
  expect_equal(result$at_skew, 1)    # (A-T)/(A+T) = (4-0)/(4+0) = 1
})


# ------------------------------------------------------------
# rb_train_classifier(): model training
# ------------------------------------------------------------

test_that("rb_train_classifier trains a model on synthetic data", {
  
  skip_if_not_installed("parsnip")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("ranger")
  skip_on_cran()
  
  set.seed(42)
  
  # Synthetic data with two distinguishable classes
  synthetic <- data.frame(
    seq_type = rep(c("12S", "COI"), each = 50),
    length = c(rnorm(50, 200, 20), rnorm(50, 600, 50)),
    gc_content = c(runif(50, 0.3, 0.5), runif(50, 0.4, 0.6)),
    at_content = c(runif(50, 0.5, 0.7), runif(50, 0.4, 0.6)),
    gc_skew = runif(100, -0.2, 0.2),
    at_skew = runif(100, -0.2, 0.2),
    stringsAsFactors = FALSE
  )
  
  model_result <- rb_train_classifier(
    synthetic,
    label_col = "seq_type",
    trees = 50,
    mtry = 3,
    min_n = 2,
    prop = 0.7,
    seed = 42
  )
  
  expect_type(model_result, "list")
  expect_true("model" %in% names(model_result))
  expect_true("metrics" %in% names(model_result))
  expect_true("test_results" %in% names(model_result))
  expect_true("label_col" %in% names(model_result))
  expect_true("features" %in% names(model_result))
})


test_that("rb_train_classifier errors on missing label column", {
  
  skip_if_not_installed("parsnip")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_on_cran()
  
  synthetic <- data.frame(
    length = rnorm(10, 500, 50),
    gc_content = runif(10),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_train_classifier(synthetic, label_col = "nonexistent_col"),
    "Missing required columns"
  )
})


# ------------------------------------------------------------
# rb_classify_sequences(): prediction
# ------------------------------------------------------------

test_that("rb_classify_sequences returns correct output structure", {
  
  skip_if_not_installed("parsnip")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("ranger")
  skip_on_cran()
  
  set.seed(42)
  
  # Train a simple model
  synthetic <- data.frame(
    seq_type = rep(c("12S", "COI"), each = 50),
    length = c(rnorm(50, 200, 20), rnorm(50, 600, 50)),
    gc_content = c(runif(50, 0.3, 0.5), runif(50, 0.4, 0.6)),
    at_content = c(runif(50, 0.5, 0.7), runif(50, 0.4, 0.6)),
    gc_skew = runif(100, -0.2, 0.2),
    at_skew = runif(100, -0.2, 0.2),
    sequence = replicate(100, paste0(sample(c("A","C","G","T"), 100, replace = TRUE), collapse = "")),
    stringsAsFactors = FALSE
  )
  
  model_result <- rb_train_classifier(
    synthetic,
    label_col = "seq_type",
    trees = 50,
    mtry = 3,
    min_n = 2,
    prop = 0.7,
    seed = 42
  )
  
  # Classify new sequences
  new_data <- data.frame(
    sequence = c(
      paste0(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = ""),
      paste0(sample(c("A","C","G","T"), 600, replace = TRUE), collapse = "")
    ),
    stringsAsFactors = FALSE
  )
  
  result <- rb_classify_sequences(model_result, new_data)
  
  expect_s3_class(result, "data.frame")
  expect_true("predicted_class" %in% names(result))
  expect_true("confidence_score" %in% names(result))
  expect_true("final_type" %in% names(result))
  expect_true("prediction_source" %in% names(result))
  expect_equal(nrow(result), 2)
})


test_that("rb_classify_sequences applies fallback for short sequences", {
  
  skip_if_not_installed("parsnip")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("ranger")
  skip_on_cran()
  
  set.seed(42)
  
  synthetic <- data.frame(
    seq_type = rep(c("12S", "COI"), each = 50),
    length = c(rnorm(50, 200, 20), rnorm(50, 600, 50)),
    gc_content = c(runif(50, 0.3, 0.5), runif(50, 0.4, 0.6)),
    at_content = c(runif(50, 0.5, 0.7), runif(50, 0.4, 0.6)),
    gc_skew = runif(100, -0.2, 0.2),
    at_skew = runif(100, -0.2, 0.2),
    sequence = replicate(100, paste0(sample(c("A","C","G","T"), 100, replace = TRUE), collapse = "")),
    stringsAsFactors = FALSE
  )
  
  model_result <- rb_train_classifier(
    synthetic,
    label_col = "seq_type",
    trees = 50,
    mtry = 3,
    min_n = 2,
    prop = 0.7,
    seed = 42
  )
  
  # Very short sequence should get fallback label
  short_data <- data.frame(
    sequence = "ACGT",  # 4 bp, below min_length = 100
    stringsAsFactors = FALSE
  )
  
  result <- rb_classify_sequences(model_result, short_data, min_length = 100)
  
  expect_equal(result$final_type, "other")
  expect_equal(result$prediction_source, "ML_low_confidence")
})


test_that("rb_classify_sequences uses fallback for low confidence predictions", {
  
  skip_if_not_installed("parsnip")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("ranger")
  skip_on_cran()
  
  set.seed(42)
  
  synthetic <- data.frame(
    seq_type = rep(c("12S", "COI"), each = 50),
    length = c(rnorm(50, 200, 20), rnorm(50, 600, 50)),
    gc_content = c(runif(50, 0.3, 0.5), runif(50, 0.4, 0.6)),
    at_content = c(runif(50, 0.5, 0.7), runif(50, 0.4, 0.6)),
    gc_skew = runif(100, -0.2, 0.2),
    at_skew = runif(100, -0.2, 0.2),
    sequence = replicate(100, paste0(sample(c("A","C","G","T"), 100, replace = TRUE), collapse = "")),
    stringsAsFactors = FALSE
  )
  
  model_result <- rb_train_classifier(
    synthetic,
    label_col = "seq_type",
    trees = 50,
    mtry = 3,
    min_n = 2,
    prop = 0.7,
    seed = 42
  )
  
  # Set confidence threshold impossibly high
  new_data <- data.frame(
    sequence = paste0(sample(c("A","C","G","T"), 300, replace = TRUE), collapse = ""),
    stringsAsFactors = FALSE
  )
  
  result <- rb_classify_sequences(
    model_result,
    new_data,
    confidence_threshold = 0.9999  # impossibly high
  )
  
  expect_equal(result$final_type, "other")
  expect_equal(result$prediction_source, "ML_low_confidence")
})


# ------------------------------------------------------------
# QC Round-trip: classifier output → QC Pass 2
# ------------------------------------------------------------

test_that("QC pipeline handles reclassified sequences correctly", {
  
  skip_if_not_installed("Biostrings")
  skip_on_cran()
  
  # Simulate sequences that were initially "other" but got reclassified
  reclassified_data <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3"),
    species = c("SpeciesA", "SpeciesA", "SpeciesB"),
    source = "ncbi",
    seq_type = c("COI", "12S", "other"),  # reclassified by ML
    sequence = c(
      strrep("ATGCGTACGTAGCTAGCTAGCTGATCGATCG", 5),   # 150bp COI
      strrep("ACGTACGTACGTACGTACGT", 8),               # 160bp 12S
      strrep("TTTTGGGGCCCC", 10)                        # 120bp still "other"
    ),
    unique_code = c("seq_00001", "seq_00002", "seq_00003"),
    is_contaminant = c(FALSE, FALSE, FALSE),
    is_numt = c(FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  
  # Run codon check on the COI sequence
  coi_result <- rb_check_codons(reclassified_data$sequence[1], gene = "COI")
  expect_type(coi_result, "list")
  expect_false(is.na(coi_result$has_stop))
  
  # Run rRNA check on the 12S sequence
  rrna_result <- rb_check_rrna_integrity(reclassified_data$sequence[2], gene = "12S")
  expect_type(rrna_result, "list")
  expect_false(is.na(rrna_result$has_short))
  expect_false(is.na(rrna_result$has_gaps))
  
  # "other" gene should still return NA for gene-specific checks
  other_codon <- rb_check_codons(reclassified_data$sequence[3], gene = "other")
  expect_true(is.na(other_codon$has_stop))
})


test_that("rb_compile_qc_flags handles mixed reclassified and original data", {
  
  # Simulate data after QC Pass 2 with mixed gene types
  mixed_data <- data.frame(
    is_contaminant = c(FALSE, FALSE, TRUE),
    is_numt = c(FALSE, TRUE, FALSE),
    has_stop = c(FALSE, NA, NA),
    frameshifted = c(FALSE, NA, NA),
    has_short = c(NA, FALSE, NA),
    has_gaps = c(NA, FALSE, NA),
    genome_flag = NA_character_,
    div_result = NA_character_,
    gene = c("COI", "12S", "other"),
    stringsAsFactors = FALSE
  )
  
  # Define extra_flag_fn for value-based flags
  extra_flags <- function(data) {
    vapply(seq_len(nrow(data)), function(i) {
      flags <- c()
      if (data$gene[i] == "other") flags <- c(flags, "manual_review_needed")
      paste(flags, collapse = "|")
    }, character(1))
  }
  
  result <- rb_compile_qc_flags(mixed_data, extra_flag_fn = extra_flags)
  
  expect_s3_class(result, "data.frame")
  expect_true("qc_flag" %in% names(result))
  
  # First row: COI, no issues -> pass
  expect_equal(result$qc_flag[1], "pass")
  
  # Second row: 12S with numt -> should have numt flag
  expect_true(grepl("numt", result$qc_flag[2]))
  
  # Third row: other with contaminant -> should have both flags
  expect_true(grepl("contaminant", result$qc_flag[3]))
  expect_true(grepl("manual_review_needed", result$qc_flag[3]))
})


test_that("QC Pass 2 applies gene-specific checks to newly classified sequences", {
  
  skip_if_not_installed("Biostrings")
  
  # Before classification: all "other"
  before_class <- data.frame(
    seq_type = c("other", "other", "other"),
    sequence = c(
      strrep("ATGCGTACGTAGCTAGCTAGCTGATCGATCG", 5),
      strrep("ACGTACGTACGTACGTACGT", 8),
      strrep("TTTTGGGGCCCC", 10)
    ),
    stringsAsFactors = FALSE
  )
  
  # After classification: some reclassified
  after_class <- before_class
  after_class$seq_type <- c("COI", "12S", "other")
  
  # Gene-specific checks should now apply to reclassified sequences
  # COI should get codon check
  expect_false(is.na(rb_check_codons(after_class$sequence[1], gene = after_class$seq_type[1])$has_stop))
  
  # 12S should get rRNA check
  expect_false(is.na(rb_check_rrna_integrity(after_class$sequence[2], gene = after_class$seq_type[2])$has_short))
  
  # "other" should still be skipped
  expect_true(is.na(rb_check_codons(after_class$sequence[3], gene = after_class$seq_type[3])$has_stop))
  expect_true(is.na(rb_check_rrna_integrity(after_class$sequence[3], gene = after_class$seq_type[3])$has_short))
})