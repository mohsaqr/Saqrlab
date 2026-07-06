# Compare Two TNA Networks

Compare two TNA network models using multiple metrics including
correlation, RMSE, and edge-level differences.

## Usage

``` r
compare_networks(
  model1,
  model2,
  metrics = c("correlation", "rmse", "edge_diff"),
  scaling = c("none", "minmax", "zscore"),
  include_self = TRUE,
  threshold = 0.05
)
```

## Arguments

- model1:

  A TNA model object (the reference/original model).

- model2:

  A TNA model object (the comparison/simulated model).

- metrics:

  Character vector. Metrics to compute. Options: "correlation" (Pearson
  correlation), "rmse" (root mean square error), "mae" (mean absolute
  error), "edge_diff" (proportion with large differences), "cosine"
  (cosine similarity), or "all" to compute everything. Default:
  c("correlation", "rmse", "edge_diff").

- scaling:

  Character. How to scale edge weights before comparison. Options:
  "none" (no scaling), "minmax" (min-max normalization to 0-1), "zscore"
  (z-score standardization). Default: "none".

- include_self:

  Logical. Whether to include self-loops in comparison. Default: TRUE.

- threshold:

  Numeric. Threshold for considering an edge as "different" in edge_diff
  metric. Default: 0.05.

## Value

A list with class "tna_comparison" containing:

- metrics: Named list of computed metrics.

- edge_comparison: Data frame comparing edges from both models.

- summary: Character summary of the comparison.

## Details

This function extracts edge weights from both models and computes
various comparison metrics. Edge weights are extracted from the
transition matrices stored in the TNA model objects.

## See also

[`compare_centralities()`](https://pak.dynasite.org/Saqrlab/reference/compare_centralities.md)
for comparing centrality profiles,
[`compare_edge_recovery()`](https://pak.dynasite.org/Saqrlab/reference/compare_edge_recovery.md)
for edge recovery metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate original data and fit model
original_data <- simulate_sequences(trans_mat, init_probs, 20, 200)
model_original <- tna::tna(original_data)

# Simulate from fitted model and fit new model
sim_data <- simulate_sequences(
  transition_matrix = model_original$weights,
  initial_probabilities = model_original$initial,
  max_seq_length = 20, num_rows = 200
)
model_sim <- tna::tna(sim_data)

# Compare the models
comparison <- compare_networks(model_original, model_sim)
print(comparison$metrics)
} # }
```
