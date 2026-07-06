# Plot Sampling Distribution for a Single Metric

Create a visualization showing the distribution of a single metric's
values across all sampling iterations.

## Usage

``` r
plot_sampling_distribution(individual_results, category, metric)
```

## Arguments

- individual_results:

  Data frame of raw metric values from iterations, as returned by
  [`run_sampling_analysis()`](https://pak.dynasite.org/Saqrlab/reference/run_sampling_analysis.md)
  or
  [`compare_tna_models()`](https://pak.dynasite.org/Saqrlab/reference/compare_network_estimation.md).

- category:

  The category of the metric to plot (e.g., "Correlations").

- metric:

  The specific metric name to plot (e.g., "Pearson").

## Value

A ggplot object showing the distribution.

## Examples

``` r
if (FALSE) { # \dontrun{
library(tna)
model <- tna(group_regulation)
results <- run_sampling_analysis(model, iterations = 100)

# Plot Pearson correlation distribution
plot_sampling_distribution(
  results$individual,
  category = "Correlations",
  metric = "Pearson"
)
} # }
```
