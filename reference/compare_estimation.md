# Compare Model Estimation Across Simulations

Run multiple simulations comparing how well different TNA model types
recover the true transition structure, using
[`tna::compare()`](http://sonsoles.me/tna/reference/compare.md) for each
comparison. Supports any model type available in the tna package.

## Usage

``` r
compare_estimation(
  models = c("tna", "ftna"),
  n_simulations = 1000,
  n_sequences = 200,
  seq_length = 25,
  n_states = 6,
  na_range = c(0, 5),
  scaling = "minmax",
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  cores = parallel::detectCores() - 1
)
```

## Arguments

- models:

  Character vector. Model types to compare. Any model function exported
  by tna package works: "tna", "ftna", "ctna", "atna", "sna", "tsn".
  Default: c("tna", "ftna").

- n_simulations:

  Integer. Number of simulations to run. Default: 1000.

- n_sequences:

  Integer. Number of sequences per simulation. Default: 200.

- seq_length:

  Integer. Maximum sequence length. Default: 25.

- n_states:

  Integer. Number of states. Default: 6.

- na_range:

  Integer vector of length 2. Range of NAs per sequence for varying
  lengths. Default: c(0, 5).

- scaling:

  Character. Scaling for tna::compare(). Default: "minmax".

- seed:

  Integer or NULL. Random seed. Default: NULL.

- verbose:

  Logical. Print progress. Default: TRUE.

- parallel:

  Logical. Use parallel processing. Default: FALSE.

- cores:

  Integer. Number of cores for parallel. Default: detectCores() - 1.

## Value

A list containing:

- comparison:

  Side-by-side comparison of key metrics across models.

- summary:

  Data frame with mean/sd of all metrics by model type.

- raw_results:

  Data frame with all simulation results.

- ranking:

  Models ranked by Pearson correlation (best to worst).

- winner:

  Model with highest Pearson correlation.

- params:

  Parameters used for the simulation.

## Details

For each simulation:

1.  Generate random transition probabilities (ground truth)

2.  Simulate sequences with optional NAs (varying lengths)

3.  Fit all specified model types

4.  Compare each to ground truth using
    [`tna::compare()`](http://sonsoles.me/tna/reference/compare.md)

5.  Collect metrics (Pearson, Spearman, Kendall, Euclidean, etc.)

Available model types (any tna package model):

- `tna`: Standard transition network analysis (probabilities)

- `ftna`: Frequency-based TNA (raw counts)

- `ctna`: Concurrent TNA

- `atna`: Absorbing TNA

- `sna`: Sequential network analysis

- `tsn`: Time-series network

- Any other model function exported by tna package

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare 2 models (default: tna vs ftna)
results <- compare_estimation(n_simulations = 100, seed = 42)

# Compare 3 models
results <- compare_estimation(
  models = c("tna", "ftna", "ctna"),
  n_simulations = 100,
  seed = 42
)

# Compare all 4 models
results <- compare_estimation(
  models = c("tna", "ftna", "ctna", "atna"),
  n_simulations = 500,
  parallel = TRUE,
  seed = 42
)

# View results
results$comparison  # Side-by-side metrics
results$ranking     # Best to worst
results$winner      # Top performer
} # }
```
