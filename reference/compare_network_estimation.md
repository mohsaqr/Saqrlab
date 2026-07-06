# Compare Network Estimation Across TNA Model Types

Perform side-by-side comparison of multiple TNA model types (tna, atna,
ctna, ftna) using comprehensive sampling analysis. Uses a statistically
correct approach: comparing models built on independent data subsets
(sample vs remaining).

## Usage

``` r
compare_network_estimation(
  data,
  model_types = c("tna", "ftna"),
  model_scaling = NULL,
  sampling_percent = 0.3,
  iterations = 100,
  seed = NULL,
  verbose = TRUE
)

compare_tna_models(
  data,
  model_types = c("tna", "ftna"),
  model_scaling = NULL,
  sampling_percent = 0.3,
  iterations = 100,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame containing sequence data.

- model_types:

  Character vector of model types to compare. Valid types are "tna",
  "atna", "ctna", "ftna". Default: c("tna", "ftna").

- model_scaling:

  A named list specifying scaling for each model type. Example: list(tna
  = "minmax", ftna = "skip"). If NULL, uses default scaling.

- sampling_percent:

  Numeric value between 0 and 1 indicating the proportion of data to
  sample. Default: 0.3.

- iterations:

  Integer specifying the number of sampling iterations for each model.
  Default: 100.

- seed:

  Integer seed for reproducible random sampling. Default: NULL.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list containing:

- individual:

  Data frame with all iteration results including model_type column.

- aggregated:

  Data frame with aggregated metrics by model type.

- ranking:

  Models ranked by mean Pearson correlation.

- winner:

  Model with highest Pearson correlation.

- params:

  Parameters used for the comparison.

## Details

The function workflow:

1.  For each model type, builds a model from the data.

2.  Runs sampling analysis comparing sample vs remaining subsets.

3.  Combines all results for comparison.

This approach is statistically correct because it compares models built
on independent data subsets, providing a true measure of estimation
stability.

## Examples

``` r
if (FALSE) { # \dontrun{
library(tna)
data(group_regulation)

# Compare tna and ftna models
results <- compare_network_estimation(group_regulation)
results$ranking
results$winner

# Compare all model types
results <- compare_network_estimation(
  group_regulation,
  model_types = c("tna", "atna", "ctna", "ftna"),
  iterations = 100
)

# Plot results
plot(results)  # Default histogram
plot(results, metric_name = "Pearson")
} # }
```
