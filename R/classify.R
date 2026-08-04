# ============================================================
# MODULE 4 — ML sequence-type classifier (replaces edna_ml_seq_type.R)
#
# DESTINATION: R/classify.R (new file). This is a distinct capability
# (model training + prediction) not represented anywhere else in what
# you've shown me, so it earns its own file rather than being folded
# into R/qc.R.
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

# --- Test ---
# rb_sequence_features(c("AATTGGCC", "GGGGCCCC"))
# expect: length = c(8, 8); gc_content = c(0.5, 1.0); at_content = c(0.5, 0.0)
#
# rb_sequence_features("")
# expect: length = 0, all content/skew values = 0 (no divide-by-zero error,
#         confirms the total_bases > 0 guard still works on an edge case)


# ------------------------------------------------------------
# rb_train_classifier()
# Generalizes Stage 3 (prepare_ml_data) + Stage 4
# (train_sequence_classifier). `label_col` and `features` are now
# parameters instead of the hardcoded formula
# seq_type ~ length + gc_content + at_content + gc_skew + at_skew —
# class labels come from whatever's actually present in data[[label_col]],
# so a dataset using a completely different marker set (or even a
# non-marker classification task built on the same feature columns)
# works unchanged.
# ------------------------------------------------------------

rb_train_classifier <- function(data, label_col = "seq_type",
                                 features = c("length", "gc_content", "at_content",
                                              "gc_skew", "at_skew"),
                                 trees = 500, mtry = 3, min_n = 5,
                                 prop = 0.8, seed = 123) {
  rb_required_columns(data, c(label_col, features))

  data[[label_col]] <- factor(data[[label_col]])

  set.seed(seed)
  split <- rsample::initial_split(data, strata = label_col, prop = prop)
  train_data <- rsample::training(split)
  test_data <- rsample::testing(split)

  model_formula <- stats::as.formula(paste(label_col, "~", paste(features, collapse = " + ")))

  recipe <- recipes::recipe(model_formula, data = train_data) %>%
    recipes::step_normalize(recipes::all_numeric_predictors()) %>%
    recipes::step_zv(recipes::all_predictors())

  rf_spec <- parsnip::rand_forest(mtry = mtry, trees = trees, min_n = min_n) %>%
    parsnip::set_engine("ranger", importance = "permutation") %>%
    parsnip::set_mode("classification")

  rf_workflow <- workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(rf_spec)

  fit <- workflows::fit(rf_workflow, train_data)

  class_preds <- stats::predict(fit, test_data)
  prob_preds <- stats::predict(fit, new_data = test_data, type = "prob")
  test_results <- dplyr::bind_cols(class_preds, prob_preds, test_data[label_col])

  metrics <- yardstick::metrics(test_results, truth = !!rlang::sym(label_col), estimate = .pred_class)
  conf_matrix <- yardstick::conf_mat(test_results, truth = !!rlang::sym(label_col), estimate = .pred_class)

  list(
    model = fit,
    metrics = metrics,
    test_results = test_results,
    conf_matrix = conf_matrix,
    label_col = label_col,          # stored so rb_classify_sequences() knows
    features = features              # what to expect without re-guessing
  )
}

# --- Test ---
# set.seed(1)
# toy <- data.frame(
#   seq_type = rep(c("markerA", "markerB", "markerC"), each = 20),  # deliberately
#                                                                     # NOT 12S/16S/COI
#   length = stats::rnorm(60, 500, 50),
#   gc_content = stats::runif(60), at_content = stats::runif(60),
#   gc_skew = stats::runif(60, -1, 1), at_skew = stats::runif(60, -1, 1)
# )
# model <- rb_train_classifier(toy, label_col = "seq_type")
# expect: trains without error; model$test_results has .pred_markerA /
#         .pred_markerB / .pred_markerC columns (proves classes aren't
#         hardcoded anywhere in this function)


# ------------------------------------------------------------
# rb_classify_sequences()
# Generalizes Stage 5 (predict_unknown_sequences_with_confidence).
# THIS IS THE FUNCTION WITH THE CRITICAL BUG — see the BUGFIX comment
# at the max_prob computation below.
# ------------------------------------------------------------

rb_classify_sequences <- function(model_result, data, confidence_threshold = 0.8,
                                   min_length = 100, sequence_col = "sequence",
                                   fallback_label = "other") {
  rb_required_columns(data, sequence_col)

  fit <- model_result$model

  # Features recomputed via rb_sequence_features() — this also removes
  # the second, slightly-differently-named copy of the feature-extraction
  # block that existed in the original predict_unknown_sequences_with_confidence()
  # (it used `seq_length` as an interim column name where the training-side
  # block used `length`; both are gone now in favor of the one function).
  features <- rb_sequence_features(data[[sequence_col]])

  pred_class <- stats::predict(fit, new_data = features, type = "class")
  pred_prob <- stats::predict(fit, new_data = features, type = "prob")

  # BUGFIX 2026-08-04 (rb_classify_sequences, generalizing
  # predict_unknown_sequences_with_confidence()):
  #
  # Original code:
  #   max_prob = pmax(.pred_12S, .pred_16S, .pred_COI)
  #
  # This hardcodes the three YZFishDB marker class names directly into
  # the pmax() call. For any dataset whose label_col contains different
  # class names (e.g. the markerA/markerB/markerC example above, or a
  # marker set with 4+ classes), the columns .pred_12S/.pred_16S/.pred_COI
  # simply would not exist in pred_prob — this would throw
  # "object '.pred_12S' not found" immediately, or worse, silently
  # succeed with wrong values if similarly-named columns happened to
  # exist from unrelated code in the same session.
  #
  # Fixed to compute the max across EVERY .pred_* column actually
  # present in pred_prob, so it works for any number or naming of
  # classes without the caller needing to know or specify them.
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

# --- Test (regression test for the pmax fix — this is the case that
#     previously would have errored or silently misbehaved) ---
# new_data <- data.frame(sequence = c("ACGTACGTACGTACGTACGT"))
# result <- rb_classify_sequences(model, new_data)
# expect (after fix): runs without error; confidence_score is the max
#   across .pred_markerA / .pred_markerB / .pred_markerC
# BEFORE fix: this call would have thrown
#   "object '.pred_12S' not found" — the markerA/B/C dataset has no
#   such column, since the original code hardcoded YZFishDB's own
#   class names into the pmax() call.
#
# Also worth a second regression test confirming a short sequence is
# always routed to fallback_label regardless of model confidence:
# rb_classify_sequences(model, data.frame(sequence = "AC"), min_length = 100)
# expect: final_type == "other" (fallback), prediction_source == "ML_low_confidence"
