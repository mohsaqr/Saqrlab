# Calculate Edge Recovery Metrics

Calculate metrics for how well edges from an original network are
recovered in a simulated/comparison network.

## Usage

``` r
compare_edge_recovery(
  original,
  simulated,
  threshold = 0.01,
  return_edges = FALSE
)

calculate_edge_recovery(...)
```

## Arguments

- original:

  A TNA model object (the original/ground truth model).

- simulated:

  A TNA model object (the simulated/recovered model).

- threshold:

  Numeric. Minimum edge weight to consider an edge as "present".
  Default: 0.01.

- return_edges:

  Logical. Whether to return detailed edge-level results. Default:
  FALSE.

- ...:

  Arguments passed to `compare_edge_recovery`.

## Value

A list containing:

- true_positives: Number of edges correctly present in both.

- false_positives: Number of edges present in simulated but not
  original.

- false_negatives: Number of edges present in original but not
  simulated.

- true_negatives: Number of edges correctly absent in both.

- precision: TP / (TP + FP).

- recall: TP / (TP + FN), also known as sensitivity.

- f1_score: Harmonic mean of precision and recall.

- accuracy: (TP + TN) / total edges.

- jaccard: TP / (TP + FP + FN), Jaccard similarity.

- edges: (Optional) Data frame with edge-level results.

## Details

This function treats edge recovery as a binary classification problem:

- True Positive: Edge present in both original and simulated.

- False Positive: Edge present in simulated but not original.

- False Negative: Edge present in original but not simulated.

- True Negative: Edge absent in both.

An edge is considered "present" if its weight exceeds the threshold.

## See also

[`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md)
for full network comparison,
[`run_bootstrap_iteration()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_iteration.md)
for bootstrap evaluation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Calculate edge recovery
recovery <- compare_edge_recovery(model_original, model_simulated)
print(sprintf("Precision: %.2f, Recall: %.2f", recovery$precision, recovery$recall))
} # }
```
