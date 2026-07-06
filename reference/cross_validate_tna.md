# Cross-Validate TNA Model Types

Perform cross-validation by testing different TNA model types on the
same data, allowing comparison of model performance.

## Usage

``` r
cross_validate_tna(
  data,
  model_types = c("relative", "frequency", "co-occurrence"),
  sampling_percent = 0.3,
  iterations = 50,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame containing sequence data.

- model_types:

  Character vector of model types to test. Default: c("relative",
  "frequency", "co-occurrence").

- sampling_percent:

  Proportion of data for sampling. Default: 0.3.

- iterations:

  Number of iterations per model type. Default: 50.

- seed:

  Random seed for reproducibility. Default: NULL.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

A list containing cross-validation results for each model type.

## Examples

``` r
if (FALSE) { # \dontrun{
data(group_regulation, package = "tna")
cv_results <- cross_validate_tna(
  group_regulation,
  model_types = c("relative", "frequency"),
  iterations = 30
)
} # }
```
