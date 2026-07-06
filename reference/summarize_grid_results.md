# Summarize Grid Simulation Results

Summarize results from
[`run_grid_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_grid_simulation.md),
filtering by parameter ranges and computing aggregated performance
metrics. Provides detailed summaries at both the setting level and
across all selected settings.

## Usage

``` r
summarize_grid_results(
  grid_results_list,
  n_sequences_range = NULL,
  seq_length_range = NULL,
  min_na_range = NULL,
  max_na_range = NULL,
  level_context = 0.05,
  print_output = TRUE,
  print_aggregated_overall = TRUE,
  print_aggregated_edges = TRUE,
  print_settings_summary = TRUE,
  num_rows_range = NULL,
  max_seq_length_range = NULL
)

analyze_grid_results(...)
```

## Arguments

- grid_results_list:

  List output from
  [`run_grid_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_grid_simulation.md).

- n_sequences_range:

  Numeric vector of length 2 (min, max) or NULL. Filter settings by
  n_sequences.

- seq_length_range:

  Numeric vector of length 2 (min, max) or NULL. Filter settings by
  seq_length.

- min_na_range:

  Numeric vector of length 2 (min, max) or NULL. Filter settings by
  min_na.

- max_na_range:

  Numeric vector of length 2 (min, max) or NULL. Filter settings by
  max_na.

- level_context:

  Numeric. Significance level used for calculations (for context in
  output). Default: 0.05.

- print_output:

  Logical. Master switch for console printing. Default: TRUE.

- print_aggregated_overall:

  Logical. Print averaged overall performance metrics across selected
  settings. Default: TRUE.

- print_aggregated_edges:

  Logical. Print aggregated edge significance summary. Default: TRUE.

- print_settings_summary:

  Logical. Print summary table for each selected setting. Default: TRUE.

- num_rows_range:

  Deprecated. Use `n_sequences_range` instead.

- max_seq_length_range:

  Deprecated. Use `seq_length_range` instead.

- ...:

  Arguments passed to `summarize_grid_results`.

## Value

A list containing:

- n_selected:

  Number of settings matching the filter criteria.

- aggregated_summary:

  List with:

  - overall_performance: Averaged metrics across settings.

  - edge_significance: Aggregated edge-level statistics.

- selected_settings_summary_df:

  Data frame with per-setting metrics calculated from total TP/TN/FP/FN
  counts.

- compiled_individual_runs:

  List with detailed run-level data:

  - all_raw_summaries: Combined bootstrap summaries.

  - all_per_edge_performance: Combined per-edge results.

  - run_level_performance_metrics: Metrics per run.

  - setting_level_summary_stats: Mean/Median/SD of run metrics.

Returns NULL (invisibly) if no settings match.

## Details

The function performs comprehensive analysis:

**Filtering**: Selects settings where all parameters fall within
specified ranges.

**Aggregation from Input Summaries**: Extracts and averages overall
performance and edge significance from the original `aggregated_summary`
in each setting.

**Run-Level Metrics**: Recomputes TP/TN/FP/FN counts from per-edge data,
then calculates Sensitivity, Specificity, FPR, FNR, Accuracy, and MCC
per run.

**Setting-Level Metrics**: Two approaches:

1.  Mean/Median/SD of run-level metrics.

2.  Metrics computed from total counts across all runs (more robust).

The `selected_settings_summary_df` uses approach \#2 for the most
accurate overall metrics per setting.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running grid simulation
grid_results <- run_grid_simulation(...)

# Analyze all results
analysis <- summarize_grid_results(grid_results)

# Filter to specific parameter ranges
analysis_filtered <- summarize_grid_results(
  grid_results,
  n_sequences_range = c(100, 200),
  seq_length_range = c(30, 50),
  max_na_range = c(0, 5)
)

# Access the detailed run-level metrics
run_metrics <- analysis$compiled_individual_runs$run_level_performance_metrics

# Access setting-level summary
setting_summary <- analysis$selected_settings_summary_df

# Old parameter names still work
analysis_filtered <- summarize_grid_results(
  grid_results,
  num_rows_range = c(100, 200),
  max_seq_length_range = c(30, 50)
)
} # }
```
