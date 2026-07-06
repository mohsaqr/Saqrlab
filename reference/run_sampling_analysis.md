# Run Sampling Analysis on TNA Models

Perform comprehensive TNA model analysis through repeated sampling and
statistical comparison. Uses a statistically correct sampling procedure
by comparing models built on independent data subsets (sample vs
remaining).

## Usage

``` r
run_sampling_analysis(
  model,
  sampling_percent = 0.3,
  iterations = 100,
  seed = NULL,
  model_scaling = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A TNA model object created with tna() or related functions.

- sampling_percent:

  Numeric value between 0 and 1 specifying the proportion of data to use
  for the sample model. Default: 0.3.

- iterations:

  Integer specifying the number of sampling iterations. Default: 100.

- seed:

  Integer seed for reproducible random sampling. Default: NULL.

- model_scaling:

  Character string to override the model's scaling: NULL (default) uses
  original model's scaling, "skip" forces no scaling, or valid TNA
  scaling options ("minmax", "max", "rank").

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list containing:

- aggregated:

  Data frame with aggregated summary statistics by metric.

- individual:

  Data frame with raw metric values from each iteration.

- params:

  Parameters used for the analysis.

## Details

This function implements a statistically correct sampling procedure:

1.  **Data Splitting**: For each iteration, splits the original data
    into:

    - Sample set (specified percentage)

    - Remaining set (complement)

2.  **Model Building**: Builds TNA models on both datasets using the
    original model's type and configurable scaling.

3.  **Model Comparison**: Compares sample model vs remaining model using
    TNA's compare() function, providing metrics across categories:

    - Correlations (Pearson, Spearman, Kendall)

    - Dissimilarities (Euclidean, Manhattan, etc.)

    - Similarities (Cosine, Jaccard, etc.)

4.  **Result Aggregation**: Collects metrics across all iterations and
    computes summary statistics (mean, sd, median, quartiles).

This approach is statistically superior to comparing sample vs original
because it compares models built on independent data subsets, providing
a true measure of sampling variability and model stability.

## Examples

``` r
if (FALSE) { # \dontrun{
library(tna)
model <- tna(group_regulation)

# Basic analysis with defaults
results <- run_sampling_analysis(model, iterations = 50)
results$aggregated

# With custom parameters
results <- run_sampling_analysis(
  model,
  sampling_percent = 0.4,
  iterations = 100,
  seed = 42
)
} # }
```
