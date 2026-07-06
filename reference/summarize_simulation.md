# Summarize Simulation Results

Compute summary statistics for simulation results, optionally grouped by
parameter combinations or other factors.

## Usage

``` r
summarize_simulation(
  results,
  by = NULL,
  metrics = c("mean", "sd", "ci"),
  value_cols = NULL,
  na.rm = TRUE
)
```

## Arguments

- results:

  A data frame of simulation results, typically from
  [`run_network_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_network_simulation.md)
  or
  [`run_bootstrap_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_simulation.md).

- by:

  Character vector. Column names to group by when computing summaries.
  Default: NULL (compute overall summary).

- metrics:

  Character vector. Which summary metrics to compute. Options: "mean"
  (arithmetic mean), "sd" (standard deviation), "median" (median value),
  "ci" (95 percent confidence interval), "min" (minimum), "max"
  (maximum), "n" (count), or "all" to compute all metrics. Default:
  c("mean", "sd", "ci").

- value_cols:

  Character vector. Names of columns containing numeric values to
  summarize. If NULL, auto-detects numeric columns. Default: NULL.

- na.rm:

  Logical. Whether to remove NA values when computing statistics.
  Default: TRUE.

## Value

A data frame with summary statistics. If `by` is specified, contains one
row per unique combination of grouping variables.

## Details

This function provides flexible summarization of simulation results. It
automatically detects numeric columns and computes requested statistics.

Confidence intervals are computed as mean +/- 1.96 \* SE, where SE = sd
/ sqrt(n).

## See also

[`summarize_networks()`](https://pak.dynasite.org/Saqrlab/reference/summarize_networks.md)
for summarizing network metrics,
[`run_network_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_network_simulation.md)
for running simulations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run a simulation
results <- run_network_simulation(trans_mat, init_probs,
                                  num_rows_values = c(50, 100, 200),
                                  n_runs = 20)

# Overall summary
summarize_simulation(results)

# Summary by sample size
summarize_simulation(results, by = "num_rows")

# Summary with all metrics
summarize_simulation(results, by = "num_rows", metrics = "all")
} # }
```
