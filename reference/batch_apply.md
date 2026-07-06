# Apply Function to Multiple Models or Datasets

Apply a function to a list of TNA models or datasets, with optional
parallel processing.

## Usage

``` r
batch_apply(
  object_list,
  fun,
  parallel = FALSE,
  cores = NULL,
  progress = TRUE,
  simplify = FALSE,
  ...
)
```

## Arguments

- object_list:

  A list of TNA model objects or data frames.

- fun:

  A function to apply to each element.

- parallel:

  Logical. Whether to use parallel processing. Default: FALSE.

- cores:

  Integer or NULL. Number of cores for parallel processing. Default:
  NULL.

- progress:

  Logical. Whether to show progress messages. Default: TRUE.

- simplify:

  Logical. Whether to simplify results if possible. Default: FALSE.

- ...:

  Additional arguments passed to `fun`.

## Value

A list of results from applying `fun` to each element. If
`simplify = TRUE` and results are atomic, returns a simplified vector.

## Details

This is a general-purpose function for batch operations on model lists.
It handles errors gracefully, returning NULL for failed operations.

## See also

[`batch_fit_models`](https://pak.dynasite.org/Saqrlab/reference/batch_fit_models.md)
for fitting multiple models.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract transition matrices from multiple models
trans_mats <- batch_apply(models, extract_transition_matrix)

# Compare each model to a reference
ref_model <- models[[1]]
comparisons <- batch_apply(
  models[-1],
  function(m) compare_networks(ref_model, m)$metrics$correlation
)

# Extract centralities with simplification
centralities <- batch_apply(
  models,
  function(m) tna::centralities(m),
  simplify = FALSE
)
} # }
```
