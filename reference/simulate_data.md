# Simulate Ready-to-Use Statistical Datasets

Generate synthetic datasets suitable for common statistical analyses.
Each dataset type contains real signal so the intended analysis works
out of the box. Different seeds produce structurally different datasets
(varying n, effect sizes, number of groups/variables).

## Usage

``` r
simulate_data(type, seed = NULL, complexity = "clean", ..., n_batch = NULL)
```

## Arguments

- type:

  Character string specifying the dataset type. One of:

  `"ttest"`

  :   Two-group comparison. Columns: `group` (factor), `score`
      (numeric).

  `"anova"`

  :   Multi-group comparison (3–5 groups). Columns: `group` (factor),
      `score` (numeric).

  `"correlation"`

  :   Correlated variables (4–7 columns). All numeric `x1`–`xp`.

  `"clusters"`

  :   Cluster structure. Numeric `x1`–`xd` plus integer `true_cluster`.

  `"factor_analysis"`

  :   Latent factor structure. All numeric `x1`–`xp`. Attributes
      `n_factors` and `loadings`.

  `"prediction"`

  :   Regression dataset. Columns: `y`, `x1`–`x4`, `cat1`, `cat2`.

  `"mlvar"`

  :   Multilevel VAR panel data. Columns: `id`, `day`, `beep`,
      `V1`–`Vd`. Attributes `true_temporal`, `true_contemporaneous`,
      `vars`.

  `"batch"`

  :   Wildcard: generates all 7 types. Returns a named list with
      `n_batch` (default 1000) datasets per type.

- seed:

  Integer or NULL. Random seed. Also determines structural parameters
  (n, effect sizes, etc.). Default: NULL.

- complexity:

  Character. Controls edge-case injection for stress testing.

  `"clean"`

  :   (default for single datasets) No edge cases.

  `"auto"`

  :   (default in batch mode) Randomly injects 0–3 edge cases per
      dataset (seed-driven, fully reproducible).

  character vector

  :   Inject specific cases. One or more of: `"na"`, `"outliers"`,
      `"ties"`, `"duplicates"`, `"constant_col"`, `"all_na_col"`,
      `"tiny_n"`, `"heavy_tailed"`, `"heteroscedastic"`,
      `"extreme_imbalance"`, `"multicollinear"`.

- ...:

  Optional overrides for structural parameters (n, n_groups, etc.).

- n_batch:

  Integer or NULL. When provided, returns a list of `n_batch` datasets
  instead of a single dataset. When `type = "batch"`, defaults to 1000L.
  Each item has `seed`, `batch_id`, and `complexity` attributes.

## Value

A `data.frame` (single dataset) or a list (when `n_batch` is non-NULL or
`type = "batch"`). Each dataset has attributes `type`, `info`, and
`complexity`. Batch items additionally have `seed` and `batch_id`.

## See also

[`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md),
[`simulate_igraph()`](https://pak.dynasite.org/Saqrlab/reference/simulate_igraph.md)

## Examples

``` r
# Single clean dataset
d <- simulate_data("ttest", seed = 42)

# With edge cases
d <- simulate_data("correlation", seed = 1, complexity = "auto")
d <- simulate_data("ttest", seed = 1, complexity = c("na", "outliers"))

# Batch of 100 ttest datasets
batch <- simulate_data("ttest", seed = 1, n_batch = 100)

# All types, 1000 datasets each
all_batches <- simulate_data("batch", seed = 1)
```
