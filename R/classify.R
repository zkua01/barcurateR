# ============================================================
# MODULE 4 — ML sequence-type classifier: feature extraction, training,
# and confidence-scored prediction. (replaces edna_ml_seq_type.R) 
# ============================================================


# ------------------------------------------------------------
# rb_sequence_features()
# Generalizes the GC/AT content+skew feature block. This logic was
# ALREADY DUPLICATED in the original script — once inside
# sequential_feature_extraction() (for known/training sequences) and
# again, nearly identically, inside the batch loop of
# predict_unknown_sequences_with_confidence() (for unknown sequences,
# with `seq_length` instead of `length` as the interim variable name).
# Pulling it into one function removes that duplication, same class of
# fix as Module 3's two-script consolidation.
#
# No manual batching/gc() loop here: str_count() is already vectorized
# over the whole input, so the original hand-rolled batch_size <- 5000
# loop is unnecessary. If you're ever working with more sequences than
# fit in memory at once, wrap calls to this function in your own
# chunking (e.g. split(sequence, ceiling(seq_along(sequence)/5000)))
# rather than baking batching into the function itself.

#' Denominator (gc_count + at_count) deliberately EXCLUDES ambiguous
#' bases (N and other IUPAC codes) — this is a training-feature
#' decision, kept as-is regardless of any separate QC-side ambiguity
#' checks (see rb_check_ambiguous_content() in R/qc.R, which uses full
#' sequence length instead, for a different purpose).
# ------------------------------------------------------------
rb_sequence_features <- function(sequence) {
  seq_upper <- toupper(sequence)

  length_val <- nchar(seq_upper)
  gc_count <- stringr::str_count(seq_upper, "[GC]")
  at_count <- stringr::str_count(seq_upper, "[AT]")
  g_count <- stringr::str_count(seq_upper, "G")
  c_count <- stringr::str_count(seq_upper, "C")
  a_count <- stringr::str_count(seq_upper, "A")
  t_count <- stringr::str_count(seq_upper, "T")

  total_bases <- gc_count + at_count
  gc_content <- ifelse(total_bases > 0, gc_count / total_bases, 0)
  at_content <- ifelse(total_bases > 0, at_count / total_bases, 0)
  gc_skew <- ifelse((g_count + c_count) > 0, (g_count - c_count) / (g_count + c_count), 0)
  at_skew <- ifelse((a_count + t_count) > 0, (a_count - t_count) / (a_count + t_count), 0)

  data.frame(
    length = length_val,
    gc_content = gc_content,
    at_content = at_content,
    gc_skew = gc_skew,
    at_skew = at_skew
  )
}



# ------------------------------------------------------------
# rb_train_classifier() : train a random forest sequence-type classifier
# Generalizes Stage 3 (prepare_ml_data) + Stage 4
# (train_sequence_classifier). `label_col` and `features` are now
# parameters instead of the hardcoded formula
# seq_type ~ length + gc_content + at_content + gc_skew + at_skew —
# class labels come from whatever's actually present in data[[label_col]],
# so a dataset using a completely different marker set (or even a
# non-marker classification task built on the same feature columns)
# works unchanged.

#'
#' @param label_col Column holding class labels — not hardcoded to any
#'   specific marker set.
#' @param features Feature columns to use, default the 5 computed by
#'   rb_sequence_features().
#'
#' Requires rsample >= 1.0.0 (relies on `strata` accepting a plain
#' string). Confirm this version against your actual target
#' environment and adjust the guard/DESCRIPTION requirement if needed.
# ------------------------------------------------------------

rb_train_classifier <- function(data,
                                label_col = "seq_type",
                                features = c(
                                  "length",
                                  "gc_content",
                                  "at_content",
                                  "gc_skew",
                                  "at_skew"
                                ),
                                trees = 500,
                                mtry = 3,
                                min_n = 5,
                                prop = 0.8,
                                seed = 123) {
  
  rb_required_columns(data, c(label_col, features))
  
  data[[label_col]] <- factor(data[[label_col]])
  
  if (nlevels(data[[label_col]]) < 2) {
    stop(
      "At least two classes are required to train a sequence classifier.",
      call. = FALSE
    )
  }
  
  set.seed(seed)
  
  split <- rsample::initial_split(data, strata = !!rlang::sym(label_col), prop = prop)
  train_data <- rsample::training(split)
  test_data <- rsample::testing(split)
  
  model_formula <- stats::as.formula(
    paste(label_col, "~", paste(features, collapse = " + "))
  )
  
  recipe <- recipes::recipe(model_formula, data = train_data) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors())
  
  rf_spec <- parsnip::rand_forest(
    mtry = mtry,
    trees = trees,
    min_n = min_n
  ) %>%
    parsnip::set_engine("ranger", importance = "permutation") %>%
    parsnip::set_mode("classification")
  
  rf_workflow <- workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(rf_spec)
  
  fit <- workflows::fit(rf_workflow, train_data)
  
  class_preds <- stats::predict(fit, test_data)
  prob_preds <- stats::predict(fit, new_data = test_data, type = "prob")
  
  test_results <- dplyr::bind_cols(
    class_preds,
    prob_preds,
    test_data[label_col]
  )
  
  metrics <- NULL
  conf_matrix <- NULL
  
  # Only compute metrics if the test set contains at least two classes.
  if (nlevels(droplevels(test_data[[label_col]])) >= 2) {
    metrics <- tryCatch(
      yardstick::metrics(
        test_results,
        truth = !!rlang::sym(label_col),
        estimate = .pred_class
      ),
      error = function(e) NULL
    )
    
    conf_matrix <- tryCatch(
      yardstick::conf_mat(
        test_results,
        truth = !!rlang::sym(label_col),
        estimate = .pred_class
      ),
      error = function(e) NULL
    )
  }
  
  list(
    model = fit,
    metrics = metrics,
    test_results = test_results,
    conf_matrix = conf_matrix,
    label_col = label_col,
    features = features
  )
}


# ------------------------------------------------------------
# rb_classify_sequences(): 
#' Classify sequences using a trained model, with confidence filtering
# Generalizes Stage 5 (predict_unknown_sequences_with_confidence).

#' @param model_result Output of rb_train_classifier().
#' @param fallback_label Label assigned when confidence is below
#'   `confidence_threshold` or the sequence is shorter than
#'   `min_length` (default "other").
# ------------------------------------------------------------
rb_classify_sequences <- function(model_result, data, confidence_threshold = 0.8,
                                   min_length = 100, sequence_col = "sequence",
                                   fallback_label = "other") {
  rb_required_columns(data, sequence_col)

  fit <- model_result$model
  features <- rb_sequence_features(data[[sequence_col]])

  pred_class <- stats::predict(fit, new_data = features, type = "class")
  pred_prob <- stats::predict(fit, new_data = features, type = "prob")

  # max_prob computed across EVERY .pred_* column actually present,
  # rather than hardcoding class names (e.g. .pred_12S/.pred_16S/
  # .pred_COI) — required for this to work on any dataset's class set.
  prob_cols <- pred_prob[, grepl("^\\.pred_", names(pred_prob)), drop = FALSE]
  max_prob <- do.call(pmax, prob_cols)

  final_label <- ifelse(
    nchar(data[[sequence_col]]) < min_length, fallback_label,
    ifelse(max_prob < confidence_threshold, fallback_label, as.character(pred_class$.pred_class))
  )

  data.frame(
    data,
    predicted_class = as.character(pred_class$.pred_class),
    confidence_score = max_prob,
    final_type = final_label,
    prediction_source = ifelse(final_label == fallback_label, "ML_low_confidence", "ML_high_confidence"),
    stringsAsFactors = FALSE
  )
}
