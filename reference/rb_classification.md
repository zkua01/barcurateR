# Machine learning classification of DNA sequences

Functions to extract sequence features, train a Random Forest
classifier, and predict marker types for unknown sequences based on
sequence composition.

- `rb_sequence_features()`: Extracts numerical features (length, GC/AT
  content, GC/AT skew) from DNA sequences.

- `rb_train_classifier()`: Trains a Random Forest model to classify
  sequences into marker types using the `tidymodels` framework.

- `rb_classify_sequences()`: Predicts marker types for new sequences
  using a trained model, applying confidence and length filtering.

## Usage

``` r
rb_sequence_features(sequence)

rb_train_classifier(
  data,
  label_col = "seq_type",
  features = c("length", "gc_content", "at_content", "gc_skew", "at_skew"),
  trees = 500,
  mtry = 3,
  min_n = 5,
  prop = 0.8,
  seed = 123
)

rb_classify_sequences(
  model_result,
  data,
  confidence_threshold = 0.8,
  min_length = 100,
  sequence_col = "sequence",
  fallback_label = "other"
)
```

## Arguments

- sequence:

  A character vector of DNA sequences.

- data:

  A data frame containing the sequences and/or features.

- label_col:

  The name of the column in `data` containing the class labels (e.g.,
  marker types like "12S", "COI").

- features:

  A character vector of column names to use as predictor features.

- trees:

  The number of trees in the Random Forest.

- mtry:

  The number of predictors randomly sampled at each split.

- min_n:

  The minimum number of data points in a terminal node.

- prop:

  The proportion of data to use for the training set (the rest is used
  for testing).

- seed:

  Random seed for reproducibility of the train/test split.

- model_result:

  The output list from `rb_train_classifier()`.

- confidence_threshold:

  Minimum prediction probability required to confidently assign a class.

- min_length:

  Minimum sequence length (in base pairs) required for classification.

- sequence_col:

  The name of the column in `data` containing the raw sequences.

- fallback_label:

  The label assigned when a sequence is too short or falls below the
  `confidence_threshold` (default `"other"`).

## Value

- `rb_sequence_features()` returns a data frame of computed features
  (`length`, `gc_content`, `at_content`, `gc_skew`, `at_skew`).

- `rb_train_classifier()` returns a list with elements: `model` (fitted
  workflow), `metrics` (yardstick metrics or NULL), `test_results` (test
  set predictions), `conf_matrix` (confusion matrix or NULL),
  `label_col`, and `features`.

- `rb_classify_sequences()` returns the input `data` augmented with
  `predicted_class`, `confidence_score`, `final_type`, and
  `prediction_source` columns.

## Details

These functions require the tidymodels packages (`parsnip`, `recipes`,
`rsample`, `yardstick`, `workflows`) and the `ranger` engine. Install
them with `install.packages('tidymodels')` and
`install.packages('ranger')`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract features
features <- rb_sequence_features(refs$sequence)
features$seq_type <- refs$seq_type

# Train classifier
model <- rb_train_classifier(features, label_col = "seq_type")
model$metrics

# Classify unknown sequences
predictions <- rb_classify_sequences(model, unknown_data)
table(predictions$prediction_source)
} # }
```
