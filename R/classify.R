if (getRversion() >= "2.15.1") {
  utils::globalVariables(".pred_class")
}

#' Machine learning classification of DNA sequences
#'
#' @description
#' Functions to extract sequence features, train a Random Forest classifier,
#' and predict marker types for unknown sequences based on sequence composition.
#'
#' * `rb_sequence_features()`: Extracts numerical features (length, GC/AT content, GC/AT skew) from DNA sequences.
#' * `rb_train_classifier()`: Trains a Random Forest model to classify sequences into marker types using the `tidymodels` framework.
#' * `rb_classify_sequences()`: Predicts marker types for new sequences using a trained model, applying confidence and length filtering.
#'
#' @param sequence A character vector of DNA sequences.
#' @param data A data frame containing the sequences and/or features.
#' @param label_col The name of the column in `data` containing the class labels (e.g., marker types like "12S", "COI").
#' @param features A character vector of column names to use as predictor features.
#' @param trees The number of trees in the Random Forest.
#' @param mtry The number of predictors randomly sampled at each split.
#' @param min_n The minimum number of data points in a terminal node.
#' @param prop The proportion of data to use for the training set (the rest is used for testing).
#' @param seed Random seed for reproducibility of the train/test split.
#' @param model_result The output list from `rb_train_classifier()`.
#' @param confidence_threshold Minimum prediction probability required to confidently assign a class.
#' @param min_length Minimum sequence length (in base pairs) required for classification.
#' @param sequence_col The name of the column in `data` containing the raw sequences.
#' @param fallback_label The label assigned when a sequence is too short or falls below the `confidence_threshold` (default `"other"`).
#'
#' @return
#' * `rb_sequence_features()` returns a data frame of computed features (`length`, `gc_content`, `at_content`, `gc_skew`, `at_skew`).
#' * `rb_train_classifier()` returns a list containing the trained model (`workflows::workflow()` fit), performance metrics, and test set results.
#' * `rb_classify_sequences()` returns the input `data` augmented with `predicted_class`, `confidence_score`, `final_type`, and `prediction_source` columns.
#'
#' @name rb_classification
#' @family machine learning
NULL

#' @rdname rb_classification
#' @export
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


#' @rdname rb_classification
#' @export
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


#' @rdname rb_classification
#' @export
rb_classify_sequences <- function(model_result, data, confidence_threshold = 0.8,
                                   min_length = 100, sequence_col = "sequence",
                                   fallback_label = "other") {
  rb_required_columns(data, sequence_col)

  fit <- model_result$model
  features <- rb_sequence_features(data[[sequence_col]])

  pred_class <- stats::predict(fit, new_data = features, type = "class")
  pred_prob <- stats::predict(fit, new_data = features, type = "prob")

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
