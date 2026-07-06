# Plot Network Estimation Comparison Results

Plot method for network_estimation objects. Default is histogram with
density overlay showing mean and median values.

## Usage

``` r
# S3 method for class 'network_estimation'
plot(x, metric_name = "Pearson", plot_type = "histogram", ...)
```

## Arguments

- x:

  A network_estimation object from compare_network_estimation().

- metric_name:

  Character. Metric to plot. Default: "Pearson".

- plot_type:

  Character. Plot type: "histogram" (default), "boxplot", "density",
  "bar", or "ridgeline".

- ...:

  Additional arguments (currently unused).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- compare_network_estimation(group_regulation)
plot(results)  # Histogram of Pearson correlation
plot(results, metric_name = "Spearman")
plot(results, plot_type = "boxplot")
} # }
```
