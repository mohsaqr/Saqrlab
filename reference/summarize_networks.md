# Summarize Multiple Networks

Compute aggregate statistics across multiple TNA network models.

## Usage

``` r
summarize_networks(
  model_list,
  include = c("density", "centrality", "edges"),
  threshold = 0.01,
  centrality_measures = c("OutStrength", "InStrength")
)
```

## Arguments

- model_list:

  A list of TNA model objects.

- include:

  Character vector. What to include in the summary. Options: "density"
  (network density), "centrality" (centrality measures), "edges" (edge
  weight statistics), or "all" for everything. Default: c("density",
  "centrality", "edges").

- threshold:

  Numeric. Minimum edge weight to consider present for density
  calculation. Default: 0.01.

- centrality_measures:

  Character vector. Which centrality measures to summarize. Default:
  c("OutStrength", "InStrength").

## Value

A list containing:

- summary_table: Data frame with per-network statistics.

- aggregate: Named list of aggregate statistics across all networks.

- n_networks: Number of networks summarized.

## Details

This function is useful for summarizing results from simulation studies
where many networks are generated. It provides both per-network and
aggregate statistics.

## See also

[`summarize_simulation()`](https://pak.dynasite.org/Saqrlab/reference/summarize_simulation.md)
for simulation result summaries,
[`batch_fit_models()`](https://pak.dynasite.org/Saqrlab/reference/batch_fit_models.md)
for fitting multiple models.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate and fit multiple networks
datasets <- lapply(1:20, function(i) {
  simulate_sequences(trans_mat, init_probs, 20, 100)
})
models <- batch_fit_models(datasets)

# Summarize networks
network_summary <- summarize_networks(models)
print(network_summary$aggregate)
} # }
```
