# Run Grid Search Over Simulation Parameters

Execute bootstrap simulations across a grid of parameter combinations.
Useful for studying how different settings affect edge recovery
performance.

## Usage

``` r
run_grid_simulation(
  Model,
  stable_transitions,
  num_runs,
  n_sequences_vec = NULL,
  seq_length_vec = NULL,
  na_range_list = list(list(min = 0, max = 0)),
  stability_prob = 0.95,
  unstable_mode = "random_jump",
  unstable_random_transition_prob = 0.5,
  unstable_perturb_noise = 0.5,
  unlikely_prob_threshold = 0.1,
  include_na = TRUE,
  consistency_range = c(0.75, 1.25),
  level = 0.05,
  num_cores = parallel::detectCores() - 1,
  num_rows_vec = NULL,
  max_seq_length_vec = NULL
)
```

## Arguments

- Model:

  A TNA model object with `weights` (transition matrix) and `inits`
  (initial probabilities).

- stable_transitions:

  List of character vectors defining ground truth stable transitions.

- num_runs:

  Integer. Number of bootstrap runs per parameter combination.

- n_sequences_vec:

  Numeric vector. Values of n_sequences to test.

- seq_length_vec:

  Numeric vector. Values of seq_length to test.

- na_range_list:

  List of lists. Each inner list has `min` and `max` elements defining
  an NA range to test.

- stability_prob:

  Numeric. Fixed stability probability. Default: 0.95.

- unstable_mode:

  Character. Fixed unstable mode. Default: "random_jump".

- unstable_random_transition_prob:

  Numeric. Fixed unstable probability. Default: 0.5.

- unstable_perturb_noise:

  Numeric. Fixed perturbation noise. Default: 0.5.

- unlikely_prob_threshold:

  Numeric. Fixed unlikely threshold. Default: 0.1.

- include_na:

  Logical. Whether to include NAs. Default: TRUE.

- consistency_range:

  Numeric vector of length 2. Bootstrap consistency range. Default:
  c(0.75, 1.25).

- level:

  Numeric. Significance level. Default: 0.05.

- num_cores:

  Integer. Number of cores for parallel processing. Default:
  detectCores() - 1.

- num_rows_vec:

  Deprecated. Use `n_sequences_vec` instead.

- max_seq_length_vec:

  Deprecated. Use `seq_length_vec` instead.

## Value

A named list where each element corresponds to a parameter combination.
Names follow the pattern `nr<n_sequences>_sl<seq_length>_na<min>-<max>`.
Each element contains:

- aggregated_summary:

  Aggregated performance and edge significance.

- individual_runs:

  Per-run details.

- successful_runs:

  Number of successful runs.

- parameters:

  The parameter values used for this combination.

## Details

The function creates a full factorial grid from:

- `n_sequences_vec` x `seq_length_vec` x `na_range_list`

For each combination, it runs
[`run_bootstrap_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_simulation.md)
and stores the results along with the parameter values used.

Progress messages are printed to track execution.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a model
trans_mat <- matrix(c(
  0.6, 0.3, 0.1,
  0.2, 0.6, 0.2,
  0.1, 0.2, 0.7
), nrow = 3, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")

Model <- list(
  weights = trans_mat,
  inits = c(A = 0.33, B = 0.34, C = 0.33)
)

stable <- list(c("A", "B"), c("B", "C"))

# Run grid search
grid_results <- run_grid_simulation(
  Model = Model,
  stable_transitions = stable,
  num_runs = 20,
  n_sequences_vec = c(50, 100, 200),
  seq_length_vec = c(20, 30, 50),
  na_range_list = list(
    list(min = 0, max = 0),
    list(min = 0, max = 5),
    list(min = 5, max = 10)
  ),
  num_cores = 4
)

# Analyze results
summarize_grid_results(grid_results)

# Old parameter names still work
grid_results <- run_grid_simulation(
  Model = Model,
  stable_transitions = stable,
  num_runs = 20,
  num_rows_vec = c(50, 100, 200),
  max_seq_length_vec = c(20, 30, 50)
)
} # }
```
