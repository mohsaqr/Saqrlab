# Inject Missing Values Under a Known Mechanism (MCAR / MAR / MNAR)

Take an existing complete data.frame and inject missing values (`NA`)
into target columns under a chosen, well-defined missingness mechanism.
This is a *transformer*, not a generator: the returned object is the
same data.frame with some cells replaced by `NA`, carrying an attached
`"missing_info"` attribute that records exactly which cells were made
missing and how. Useful for testing imputation methods and
missing-data-aware estimators against a known ground truth.

The three mechanisms follow the standard Rubin taxonomy:

- `"MCAR"`:

  Missing Completely At Random. Each cell of the target columns is
  independently set missing with probability `prop`.

- `"MAR"`:

  Missing At Random. The probability that a target cell is missing
  depends on an *observed* `predictor` column, never on the
  to-be-missing value itself. Probabilities are calibrated (via a
  rank/logistic weighting on the predictor) so the average missing
  fraction is approximately `prop`.

- `"MNAR"`:

  Missing Not At Random. The probability that a target cell is missing
  depends on that cell's *own* (numeric) value: larger values are more
  likely to go missing. Calibrated so the average missing fraction per
  column is approximately `prop`.

## Usage

``` r
inject_missingness(
  data,
  mechanism = c("MCAR", "MAR", "MNAR"),
  prop = 0.1,
  cols = NULL,
  predictor = NULL,
  seed = NULL
)
```

## Arguments

- data:

  A data.frame, or a
  [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
  object (its `$data` is used). The complete data into which to inject
  missingness.

- mechanism:

  Character. One of `"MCAR"`, `"MAR"`, `"MNAR"` (matched via
  [`match.arg`](https://rdrr.io/r/base/match.arg.html)).

- prop:

  Numeric scalar in `[0, 1]`. Target average fraction of cells in the
  target columns to set missing. Default: `0.1`.

- cols:

  Character vector of target column names, or `NULL` (the default)
  meaning all columns. For `"MAR"`, the `predictor` column is
  automatically excluded from the targets.

- predictor:

  Character scalar naming the observed column that drives missingness
  under `"MAR"`. Required for `"MAR"`, ignored otherwise.

- seed:

  Integer or `NULL`. If non-`NULL`,
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called for
  reproducibility.

## Value

The input data.frame with `NA`s injected into the target columns. The
returned object additionally carries an attribute `"missing_info"`, a
list with elements:

- `mechanism`:

  The resolved mechanism string.

- `prop`:

  The requested target proportion.

- `cols`:

  The target column names actually used.

- `n_missing`:

  Total number of cells set missing.

- `realized_prop`:

  Realized fraction of target cells set missing
  (`n_missing / (nrow * length(cols))`).

- `indicator`:

  A logical matrix (`nrow(data)` x `length(cols)`) marking the cells
  that were injected as `TRUE`, with column names matching `cols`.

## Examples

``` r
df <- simulate_prediction(
  n = 500,
  coefs = c("(Intercept)" = 1, x1 = 2, x2 = -1),
  seed = 1
)$data

# MCAR: 20% of every column missing at random
mcar <- inject_missingness(df, mechanism = "MCAR", prop = 0.2, seed = 7)
attr(mcar, "missing_info")$realized_prop
#> [1] 0.2053333
str(attr(mcar, "missing_info"))
#> List of 6
#>  $ mechanism    : chr "MCAR"
#>  $ prop         : num 0.2
#>  $ cols         : chr [1:3] "y" "x1" "x2"
#>  $ n_missing    : int 308
#>  $ realized_prop: num 0.205
#>  $ indicator    : logi [1:500, 1:3] FALSE FALSE TRUE TRUE FALSE FALSE ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "y" "x1" "x2"

# MAR: missingness in y driven by observed x1
mar <- inject_missingness(df, mechanism = "MAR", prop = 0.2,
                          cols = "y", predictor = "x1", seed = 7)
colSums(is.na(mar))
#>  y x1 x2 
#> 96  0  0 

# MNAR: large y-values more likely missing
mnar <- inject_missingness(df, mechanism = "MNAR", prop = 0.2,
                           cols = "y", seed = 7)
attr(mnar, "missing_info")$n_missing
#> [1] 92
```
