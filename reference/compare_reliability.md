# Compare Reliability Across Data Conditions

Run multiple simulations assessing split-half reliability under varying
data conditions (number of sequences, sequence length, number of
states). Uses
[`tna::reliability()`](http://sonsoles.me/tna/reference/reliability.md)
for each simulation with randomly sampled parameters from user-specified
ranges.

## Usage

``` r
compare_reliability(
  n_simulations = 1000,
  n_sequences = c(50, 500),
  seq_length = c(10, 50),
  n_states = 6,
  alpha = 1,
  diag_c = 0,
  na_range = c(0, 5),
  model_type = "tna",
  reliability_iter = 100,
  reliability_split = 0.5,
  scaling = "none",
  seed = NULL,
  verbose = TRUE,
  parallel = FALSE,
  cores = parallel::detectCores() - 1
)
```

## Arguments

- n_simulations:

  Integer. Number of simulations to run. Default: 1000.

- n_sequences:

  Integer or integer vector of length 2. Number of sequences per
  simulation. If length 2, each simulation samples a random integer in
  that range. Default: c(50, 500).

- seq_length:

  Integer or integer vector of length 2. Sequence length. If length 2,
  each simulation samples a random integer in that range. Default: c(10,
  50).

- n_states:

  Integer or integer vector of length 2. Number of states. If length 2,
  each simulation samples a random integer in that range. Default: 6.

- alpha:

  Numeric or numeric vector of length 2. Dirichlet concentration
  parameter for
  [`generate_probabilities()`](https://pak.dynasite.org/Saqrlab/reference/generate_probabilities.md).
  If length 2, each simulation samples a random value from
  `runif(1, alpha[1], alpha[2])`. Small values (e.g., 0.1) produce
  sparse transition matrices; large values (e.g., 10) produce
  near-uniform matrices. Default: 1.

- diag_c:

  Numeric or numeric vector of length 2. Diagonal boost for
  [`generate_probabilities()`](https://pak.dynasite.org/Saqrlab/reference/generate_probabilities.md).
  If length 2, each simulation samples a random value from
  `runif(1, diag_c[1], diag_c[2])`. Higher values create "sticky" states
  with strong self-transitions. Default: 0.

- na_range:

  Integer vector of length 2. Range of NAs per sequence. Default: c(0,
  5).

- model_type:

  Character. Which model to fit. Default: "tna".

- reliability_iter:

  Integer. Number of iterations for
  [`tna::reliability()`](http://sonsoles.me/tna/reference/reliability.md).
  Default: 100.

- reliability_split:

  Numeric. Split proportion for
  [`tna::reliability()`](http://sonsoles.me/tna/reference/reliability.md).
  Default: 0.5.

- scaling:

  Character. Scaling for
  [`tna::reliability()`](http://sonsoles.me/tna/reference/reliability.md).
  Default: "none".

- seed:

  Integer or NULL. Random seed. Default: NULL.

- verbose:

  Logical. Print progress. Default: TRUE.

- parallel:

  Logical. Use parallel processing. Default: FALSE.

- cores:

  Integer. Number of cores for parallel. Default: detectCores() - 1.

## Value

A list of class `"tna_reliability_comparison"` containing:

- raw_results:

  Data frame with all simulation results (one row per metric per
  simulation, with parameter columns).

- summary:

  Data frame with mean/sd aggregated across simulations grouped by
  metric.

- params:

  List of input parameters for reference.

## Details

For each simulation:

1.  Randomly sample `n_sequences`, `seq_length`, `n_states` from
    user-specified ranges.

2.  Generate ground truth via
    [`generate_probabilities()`](https://pak.dynasite.org/Saqrlab/reference/generate_probabilities.md).

3.  Simulate sequences via
    [`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
    with those parameters and NAs.

4.  Fit a TNA model (e.g., `tna::tna(sequences)`).

5.  Run
    [`tna::reliability()`](http://sonsoles.me/tna/reference/reliability.md)
    on the fitted model.

6.  Extract the `$summary` table (22 metrics with
    mean/sd/median/min/max/q25/q75).

7.  Append simulation parameters as columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Quick test
res <- compare_reliability(
  n_simulations = 10,
  reliability_iter = 20,
  seed = 42
)
print(res)

# Vary number of states too
res <- compare_reliability(
  n_simulations = 100,
  n_states = c(3, 8),
  reliability_iter = 50,
  seed = 42
)
plot(res)

# Parallel execution
res <- compare_reliability(
  n_simulations = 500,
  parallel = TRUE,
  seed = 42
)
} # }
```
