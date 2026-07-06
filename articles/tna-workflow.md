# End-to-End TNA Workflow

``` r

library(Saqrlab)
```

## Overview

This vignette demonstrates complete workflows for Temporal Network
Analysis using Saqrlab and the tna package. We cover:

1.  Basic TNA workflow
2.  Model comparison
3.  Network comparison metrics
4.  Batch processing
5.  Grid simulations

## Workflow 1: Basic TNA Analysis

### Step 1: Simulate Data

``` r

# Generate sequences with metacognitive learning states
sequences <- simulate_sequences(
  n_sequences = 200,
  seq_length = 25,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

head(sequences)
#>         V1       V2       V3       V4       V5       V6       V7       V8
#> 1     Plan   Reason    Judge  Process Memorize  Process     Plan   Reason
#> 2 Retrieve Retrieve Memorize  Process     Plan  Process  Process  Process
#> 3     Plan     Plan Memorize   Reason     Plan    Judge   Reason    Judge
#> 4 Memorize    Judge   Reason    Judge   Reason     Plan     Plan  Process
#> 5 Retrieve Memorize Retrieve Memorize Retrieve Memorize Retrieve Memorize
#> 6  Process  Process    Judge   Reason  Process  Process    Judge   Reason
#>         V9     V10      V11      V12      V13      V14      V15      V16
#> 1 Retrieve   Judge Retrieve Memorize   Reason   Reason    Judge   Reason
#> 2     Plan  Reason   Reason   Reason    Judge   Reason    Judge   Reason
#> 3  Process    Plan    Judge Retrieve Memorize Retrieve Memorize   Reason
#> 4  Process   Judge  Process  Process  Process    Judge Retrieve    Judge
#> 5 Memorize  Reason    Judge Retrieve    Judge   Reason    Judge   Reason
#> 6   Reason Process  Process  Process Memorize  Process  Process Memorize
#>        V17     V18      V19      V20      V21    V22     V23      V24      V25
#> 1    Judge Process     Plan Memorize Retrieve Reason   Judge   Reason  Process
#> 2    Judge  Reason    Judge  Process    Judge Reason   Judge   Reason    Judge
#> 3    Judge  Reason    Judge   Reason   Reason  Judge Process     Plan    Judge
#> 4   Reason  Reason Retrieve    Judge   Reason Reason    Plan  Process     Plan
#> 5    Judge  Reason   Reason    Judge   Reason  Judge  Reason     Plan  Process
#> 6 Retrieve  Reason Retrieve Memorize   Reason  Judge  Reason Retrieve Memorize
dim(sequences)
#> [1] 200  25
```

### Step 2: Fit TNA Model

``` r

library(tna)

# Fit standard TNA model
model <- fit_network_model(sequences, "tna")

# View model summary
print(model)
```

### Step 3: Extract Components

``` r

# Transition matrix
trans_mat <- extract_transition_matrix(model)
print(round(trans_mat, 3))

# Initial probabilities
init_probs <- extract_initial_probs(model)
print(round(init_probs, 3))

# Edge list (filtered by threshold)
edges <- extract_edges(model, threshold = 0.05)
head(edges)
```

### Step 4: Visualize

``` r

# Plot the network
plot(model)

# Get centrality measures
centralities(model)
```

## Workflow 2: Comparing Model Types

Saqrlab supports multiple TNA model variants:

| Model  | Description                        |
|--------|------------------------------------|
| `tna`  | Standard Temporal Network Analysis |
| `ftna` | Filtered TNA (removes weak edges)  |
| `ctna` | Conditional TNA (time-varying)     |
| `atna` | Aggregated TNA                     |

### Fit Multiple Models

``` r

# Generate sequences
sequences <- simulate_sequences(
  n_sequences = 200,
  seq_length = 30,
  n_states = 6,
  seed = 42
)

# Fit all model types
models <- list(
  tna = fit_network_model(sequences, "tna"),
  ftna = fit_network_model(sequences, "ftna"),
  ctna = fit_network_model(sequences, "ctna"),
  atna = fit_network_model(sequences, "atna")
)
```

### Compare Models

``` r

# Compare each model to TNA as reference
ref_model <- models$tna

comparisons <- lapply(names(models)[-1], function(m) {
  comp <- compare_networks(ref_model, models[[m]])
  data.frame(
    model = m,
    correlation = comp$metrics$correlation,
    rmse = comp$metrics$rmse,
    edge_diff = comp$metrics$edge_diff
  )
})

do.call(rbind, comparisons)
```

## Workflow 3: Network Comparison

### Using `compare_networks()`

Compare two networks with multiple metrics:

``` r

# Generate two different networks
net1 <- simulate_tna_networks(1, n_states = 6, n_sequences = 150, seed = 42)
net2 <- simulate_tna_networks(1, n_states = 6, n_sequences = 150, seed = 123)

# Compare
comparison <- compare_networks(
  model1 = net1$network_1$model,
  model2 = net2$network_1$model,
  metrics = c("correlation", "rmse", "mae", "edge_diff", "cosine"),
  scaling = "none",
  include_self = TRUE,
  threshold = 0.05
)

# View metrics
comparison$metrics

# View edge-level comparison
head(comparison$edge_comparison)

# Print summary
cat(comparison$summary)
```

### Using `compare_centralities()`

Compare centrality profiles between networks:

``` r

cent_comp <- compare_centralities(
  model1 = net1$network_1$model,
  model2 = net2$network_1$model,
  measures = c("OutStrength", "InStrength", "Betweenness"),
  method = "both"
)

cent_comp$correlations
```

### Using `compare_edge_recovery()`

Evaluate edge recovery performance:

``` r

recovery <- compare_edge_recovery(
  original = net1$network_1$model,
  simulated = net2$network_1$model,
  threshold = 0.05,
  return_edges = TRUE
)

# Classification metrics
cat("Precision:", recovery$precision, "\n")
cat("Recall:", recovery$recall, "\n")
cat("F1 Score:", recovery$f1_score, "\n")
cat("Accuracy:", recovery$accuracy, "\n")

# Confusion matrix
cat("\nConfusion Matrix:\n")
cat("TP:", recovery$true_positives, "\n")
cat("FP:", recovery$false_positives, "\n")
cat("FN:", recovery$false_negatives, "\n")
cat("TN:", recovery$true_negatives, "\n")
```

## Workflow 4: Batch Processing

### Fitting Multiple Models

Process multiple datasets efficiently:

``` r

# Generate 20 datasets
datasets <- lapply(1:20, function(i) {
  simulate_sequences(
    n_sequences = 100,
    seq_length = 20,
    n_states = 5,
    seed = i
  )
})

# Fit models in parallel
models <- batch_fit_models(
  data_list = datasets,
  model_type = "tna",
  parallel = TRUE,
  cores = 4,
  progress = TRUE
)

length(models)
```

### Applying Functions to Multiple Models

``` r

# Extract all transition matrices
trans_matrices <- batch_apply(
  object_list = models,
  fun = extract_transition_matrix,
  parallel = TRUE
)

# Get network densities
densities <- batch_apply(
  models,
  function(m) {
    mat <- extract_transition_matrix(m)
    sum(mat > 0.05) / length(mat)
  },
  simplify = TRUE
)

summary(densities)
```

### Comparing to Reference

``` r

# Compare all models to first as reference
ref_model <- models[[1]]

correlations <- batch_apply(
  models[-1],
  function(m) compare_networks(ref_model, m)$metrics$correlation,
  simplify = TRUE
)

summary(correlations)
```

## Workflow 5: Grid Simulations

### Create Parameter Grid

``` r

# Define parameter combinations
param_grid <- generate_param_grid(
  param_ranges = list(
    n_sequences = c(50, 200),
    seq_length = c(15, 40),
    n_states = c(4, 8)
  ),
  n = 30,
  method = "lhs"
)

head(param_grid)
nrow(param_grid)
```

### Run Grid Simulation

``` r

# Run simulations across the grid
grid_results <- run_grid_simulation(
  param_grid = param_grid,
  n_runs_per_setting = 10,
  model_type = "tna",
  parallel = TRUE,
  seed = 42
)
```

### Analyze Results

``` r

# Analyze grid results
analysis <- summarize_grid_results(grid_results)

# Summary by parameter
analysis$by_n_sequences
analysis$by_seq_length
analysis$by_n_states
```

## Workflow 6: Summarizing Results

### Simulation Summary

``` r

# After running simulations, summarize results
summary_all <- summarize_simulation(
  results = grid_results,
  metrics = c("mean", "sd", "ci")
)

# Summary by specific parameter
summary_by_n <- summarize_simulation(
  results = grid_results,
  by = "n_sequences",
  metrics = "all"
)

print(summary_by_n)
```

### Network Summary

``` r

# Summarize multiple network models
network_summary <- summarize_networks(
  model_list = models,
  include = c("density", "centrality", "edges"),
  threshold = 0.01,
  centrality_measures = c("OutStrength", "InStrength")
)

# Per-network statistics
network_summary$summary_table

# Aggregate statistics
network_summary$aggregate
```

## Data Format Conversions

### Wide to Long

``` r

# Convert sequences to long format
long_data <- wide_to_long(
  data = sequences,
  id_col = NULL,           # Auto-generate IDs
  time_prefix = "V",       # Column prefix
  action_col = "Action",
  time_col = "Time",
  drop_na = TRUE
)

head(long_data)
```

### Long to Wide

``` r

# Convert back to wide format
wide_data <- long_to_wide(
  data = long_data,
  id_col = "id",
  time_col = "Time",
  action_col = "Action",
  time_prefix = "V",
  fill_na = TRUE
)

head(wide_data)
```

### Prepare for TNA Package

``` r

# Auto-detect format and prepare
tna_ready <- prepare_for_tna(
  data = my_data,
  type = "auto",
  state_names = c("Plan", "Execute", "Review"),
  validate = TRUE
)
```

## Complete Example: Recovery Study

This example simulates a complete study of network recovery across
sample sizes:

``` r

library(Saqrlab)
library(tna)

# 1. Create ground truth network
true_probs <- generate_probabilities(n_states = 5, seed = 42)

# 2. Test different sample sizes
sample_sizes <- c(25, 50, 100, 200, 500)
results <- list()

for (n in sample_sizes) {
  # Run 20 replications
  replications <- lapply(1:20, function(rep) {
    # Generate sequences from true parameters
    seqs <- simulate_sequences(
      trans_matrix = true_probs$transition_matrix,
      init_probs = true_probs$initial_probs,
      n_sequences = n,
      seq_length = 25,
      seed = rep * 1000 + n
    )

    # Fit model
    model <- fit_network_model(seqs, "tna")

    # Compare to true matrix
    estimated_mat <- extract_transition_matrix(model, type = "scaled")
    cor(as.vector(true_probs$transition_matrix), as.vector(estimated_mat))
  })

  results[[as.character(n)]] <- unlist(replications)
}

# 3. Summarize
summary_df <- data.frame(
  sample_size = sample_sizes,
  mean_correlation = sapply(results, mean),
  sd_correlation = sapply(results, sd)
)

print(summary_df)
```

## Tips for Effective Workflows

1.  **Start with known parameters**: Use
    [`generate_probabilities()`](https://pak.dynasite.org/Saqrlab/reference/generate_probabilities.md)
    to create ground truth
2.  **Use seeds consistently**: For reproducible simulation studies
3.  **Leverage parallel processing**: Use `parallel = TRUE` in batch
    functions
4.  **Build incrementally**: Test with small samples before scaling up
5.  **Track all parameters**: Use `include_params = TRUE` when
    simulating

## Summary

| Workflow | Key Functions |
|----|----|
| Basic TNA | [`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md), [`fit_network_model()`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md), `extract_*()` |
| Model comparison | [`fit_network_model()`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md) with different types |
| Network metrics | [`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md), [`compare_centralities()`](https://pak.dynasite.org/Saqrlab/reference/compare_centralities.md) |
| Batch processing | [`batch_fit_models()`](https://pak.dynasite.org/Saqrlab/reference/batch_fit_models.md), [`batch_apply()`](https://pak.dynasite.org/Saqrlab/reference/batch_apply.md) |
| Grid simulation | [`generate_param_grid()`](https://pak.dynasite.org/Saqrlab/reference/generate_param_grid.md), [`run_grid_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_grid_simulation.md) |
| Summarizing | [`summarize_simulation()`](https://pak.dynasite.org/Saqrlab/reference/summarize_simulation.md), [`summarize_networks()`](https://pak.dynasite.org/Saqrlab/reference/summarize_networks.md) |
