# Compare Centrality Profiles

Compare centrality measures between two TNA networks using correlation
or rank-based methods.

## Usage

``` r
compare_centralities(
  model1,
  model2,
  measures = c("OutStrength", "InStrength", "Betweenness"),
  method = c("both", "correlation", "rank")
)
```

## Arguments

- model1:

  A TNA model object (the reference/original model).

- model2:

  A TNA model object (the comparison/simulated model).

- measures:

  Character vector. Centrality measures to compare. Options:
  "OutStrength", "InStrength", "Betweenness", "Closeness", or "all".
  Default: c("OutStrength", "InStrength", "Betweenness").

- method:

  Character. Comparison method: "correlation" (Pearson), "rank"
  (Spearman), or "both". Default: "both".

## Value

A list containing:

- correlations: Named list of correlation values by measure.

- centrality_comparison: Data frame with centrality values for each
  state.

- summary: Character summary of the comparison.

## Details

Centrality measures are extracted from the TNA model objects using
[`tna::centralities()`](http://sonsoles.me/tna/reference/centralities.md).
The function compares how well the centrality profiles are preserved
between the original and simulated/comparison model.

## See also

[`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md)
for full network comparison.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare centrality profiles between models
comparison <- compare_centralities(model_original, model_simulated)
print(comparison$correlations)
} # }
```
