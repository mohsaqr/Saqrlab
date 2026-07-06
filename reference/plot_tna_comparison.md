# Plot TNA Model Comparison Results

Visualize the distribution of metrics from TNA model comparisons.
Supports multiple plot types including histogram (default), boxplot,
ridgeline, bar, and density plots.

## Usage

``` r
plot_tna_comparison(
  results,
  plot_type = "histogram",
  metric_name = "Pearson",
  model_types = NULL,
  binwidth = NULL
)
```

## Arguments

- results:

  A data frame containing comparison results from
  [`compare_network_estimation()`](https://pak.dynasite.org/Saqrlab/reference/compare_network_estimation.md)
  (the individual element) or
  [`run_sampling_analysis()`](https://pak.dynasite.org/Saqrlab/reference/run_sampling_analysis.md).

- plot_type:

  Character string specifying the plot type: "histogram" (default),
  "boxplot", "ridgeline", "bar", or "density".

- metric_name:

  Character string specifying a single metric to plot. Default:
  "Pearson". If NULL, all metrics will be plotted.

- model_types:

  Character vector specifying the models to include. If NULL (default),
  all models in the data will be used.

- binwidth:

  Numeric value for histogram bin width. If NULL (default), calculated
  automatically.

## Value

A ggplot object.

## Details

Available plot types:

- boxplot:

  Box plots showing distribution by model type.

- ridgeline:

  Ridge line plots for density comparison.

- bar:

  Bar plots with mean values and error bars.

- histogram:

  Histogram with density overlay (requires metric_name).

- density:

  Density plots with mean/median lines (requires metric_name).

## Examples

``` r
if (FALSE) { # \dontrun{
library(tna)
data(group_regulation)

# Run comparison
results <- compare_tna_models(
  group_regulation,
  model_types = c("tna", "atna"),
  iterations = 30
)

# Boxplot of Pearson correlation
plot_tna_comparison(results$individual, "boxplot", metric_name = "Pearson")

# Histogram of single metric
plot_tna_comparison(results$individual, "histogram", metric_name = "Pearson")

# All metrics as bar chart
plot_tna_comparison(results$individual, "bar")
} # }
```
