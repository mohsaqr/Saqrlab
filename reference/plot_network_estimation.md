# Plot Network Estimation Results

Create publication-ready plots for network estimation comparison
results. Default is a histogram optimized for high iteration counts.

## Usage

``` r
plot_network_estimation(
  results,
  metric_name = "Pearson",
  plot_type = "histogram",
  show_stats = TRUE,
  bins = NULL
)
```

## Arguments

- results:

  A network_estimation object or data frame with individual results.

- metric_name:

  Character. Metric to plot. Default: "Pearson".

- plot_type:

  Character. One of "histogram" (default), "boxplot", "density", "bar",
  or "ridgeline".

- show_stats:

  Logical. Show mean/median statistics. Default: TRUE.

- bins:

  Integer. Number of bins for histogram. Default: auto-calculated based
  on iterations (uses Sturges' rule).

## Value

A ggplot object.

## Details

The histogram plot is optimized for high iteration counts (100+) and
shows:

- Turquoise bars with density overlay

- Green vertical line for mean

- Orange dotted line for median

- Subtitle with exact mean and median values

## Examples

``` r
if (FALSE) { # \dontrun{
results <- compare_network_estimation(group_regulation, iterations = 100)

# Default histogram
plot_network_estimation(results)

# Different metrics
plot_network_estimation(results, metric_name = "Spearman")
plot_network_estimation(results, metric_name = "Euclidean")

# Different plot types
plot_network_estimation(results, plot_type = "boxplot")
plot_network_estimation(results, plot_type = "density")
} # }
```
