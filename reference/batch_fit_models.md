# Fit Models to Multiple Datasets

Fit TNA models to multiple datasets at once, with optional parallel
processing.

## Usage

``` r
batch_fit_models(
  data_list,
  model_type = c("tna", "ftna", "ctna", "atna"),
  parallel = FALSE,
  cores = NULL,
  progress = TRUE,
  ...
)
```

## Arguments

- data_list:

  A list of data frames, each containing sequence data.

- model_type:

  Character. Type of model to fit: "tna", "ftna", "ctna", or "atna".
  Default: "tna".

- parallel:

  Logical. Whether to use parallel processing. Default: FALSE.

- cores:

  Integer or NULL. Number of cores to use for parallel processing. If
  NULL, uses `parallel::detectCores() - 1`. Ignored if
  `parallel = FALSE`. Default: NULL.

- progress:

  Logical. Whether to show progress messages. Default: TRUE.

- ...:

  Additional arguments passed to the model fitting function.

## Value

A list of fitted model objects, with the same names/indices as
`data_list`. Failed fits are returned as NULL with a warning.

## Details

This function provides a convenient way to fit TNA models to many
datasets, such as when running simulations or analyzing multiple groups.

For parallel processing on Windows, set up a parallel backend first
using the `future` package.

## See also

[`fit_network_model`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md)
for fitting a single model,
[`batch_apply`](https://pak.dynasite.org/Saqrlab/reference/batch_apply.md)
for applying functions to model lists.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate multiple datasets
datasets <- lapply(1:10, function(i) {
  simulate_sequences(trans_mat, init_probs, max_seq_length = 20, num_rows = 100)
})

# Fit models in sequence
models <- batch_fit_models(datasets, model_type = "tna")

# Fit models in parallel
models <- batch_fit_models(datasets, model_type = "tna",
                           parallel = TRUE, cores = 4)

# Check results
sapply(models, function(m) if (!is.null(m)) "OK" else "Failed")
} # }
```
