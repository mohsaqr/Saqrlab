# Simulate Longitudinal Panel Data

Generate synthetic multilevel time-series data suitable for testing
longitudinal models: multilevel VAR (mlVAR), random-intercept
cross-lagged panel models (RI-CLPM), latent growth curves, and similar.

The data-generating process is a VAR(1) model with three layers:

- Temporal:

  Autoregressive and cross-lagged effects via matrix `B`. At each time
  point, `y(t) = mu_i + B (y(t-1) - mu_i) + innovation`.

- Contemporaneous:

  Within-time correlations among innovations, given by the Cholesky of
  `contemporaneous`.

- Between-person:

  Person-specific means `mu_i` drawn from `N(grand_means, between)`.

## Usage

``` r
simulate_longitudinal(
  n = 50L,
  tp = 50L,
  vars = 4L,
  temporal = NULL,
  contemporaneous = NULL,
  between = NULL,
  grand_means = NULL,
  innovation_sd = 1,
  beeps_per_day = NULL,
  ar_range = c(0.2, 0.5),
  cross_range = c(0.05, 0.2),
  n_cross = NULL,
  complexity = "clean",
  seed = NULL
)
```

## Arguments

- n:

  Integer. Number of subjects (persons). Default 50.

- tp:

  Integer. Number of time points per subject. Default 50.

- vars:

  Integer or character vector. If integer, the number of variables
  (named `V1`...`Vp`). If character, used as variable names directly.
  Default 4.

- temporal:

  Numeric matrix (`p x p`). The VAR(1) coefficient matrix `B`. Element
  `[i,j]` is the effect of variable `j` at `t-1` on variable `i` at `t`.
  Eigenvalues should have modulus \< 1 for stationarity. Default:
  diagonal 0.3 (autoregressive only, no cross-lags).

- contemporaneous:

  Numeric matrix (`p x p`). Correlation matrix of the innovation terms
  (within-time, within-person). Must be symmetric and positive-definite.
  Default: identity (uncorrelated innovations).

- between:

  Numeric matrix (`p x p`). Covariance matrix of the person-specific
  means. Controls between-person variability. Default: identity.

- grand_means:

  Numeric vector of length `p`. Population-level means for each
  variable. Default: all zeros.

- innovation_sd:

  Numeric scalar or vector of length `p`. Standard deviation(s) of the
  innovation terms (before applying contemporaneous correlation).
  Default: 1.

- beeps_per_day:

  Integer or NULL. If non-NULL, adds `day` and `beep` columns (ESM/EMA
  structure). `tp` must be divisible by `beeps_per_day`. Default: NULL
  (single `time` column).

- ar_range:

  Numeric vector of length 2. When `temporal` is NULL, diagonal
  (autoregressive) coefficients are drawn uniformly from this range.
  Default: `c(0.2, 0.5)`.

- cross_range:

  Numeric vector of length 2. When `temporal` is NULL, off-diagonal
  (cross-lag) coefficients are drawn uniformly from this range (with
  random sign). Default: `c(0.05, 0.2)`.

- n_cross:

  Integer. When `temporal` is NULL, the number of non-zero cross-lagged
  effects to include. Default: `min(3, p*(p-1))`.

- complexity:

  Character. Edge-case injection, same semantics as
  [`simulate_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md).
  One of `"clean"` (default), `"auto"`, or a character vector of
  specific cases (e.g. `c("na", "outliers")`). Supported cases: `"na"`,
  `"outliers"`, `"heavy_tailed"`, `"heteroscedastic"`, `"tiny_n"`.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `id`, `time` (or `day`/`beep`), and
  `V1`...`Vp`.

- `$params`:

  list with `temporal` (B matrix), `contemporaneous` (innovation
  correlation), `between` (person-mean covariance), `grand_means`,
  `innovation_sd`, `n`, `tp`, `var_names`.

## Details

The temporal matrix `B` must be stationary (all eigenvalues inside the
unit circle). If an auto-generated or user-supplied `B` violates this,
it is rescaled by `0.95 / max(Mod(eigen(B)$values))`.

For ESM-style data with `beeps_per_day`, the first observation of each
day is treated as a new starting point (no carry-over from previous
day's last beep), matching the assumption in most ESM-VAR software.

## Examples

``` r
# Basic: 30 subjects, 60 time points, 3 variables
r <- simulate_longitudinal(n = 30, tp = 60, vars = 3, seed = 42)
head(r$data)
#>   id time        V1         V2         V3
#> 1  1    1 0.3797833 -2.2436730  0.4158307
#> 2  1    2 1.6277822 -0.4671853  1.2661451
#> 3  1    3 0.5634720  1.2353289 -0.0260383
#> 4  1    4 1.1672413 -0.3871001  0.4766459
#> 5  1    5 1.5363739 -0.5293727  0.4341972
#> 6  1    6 1.6313855 -1.1431780 -0.1041086
r$params$temporal  # the true VAR(1) matrix
#>             V1        V2        V3
#> V1  0.47444181 0.0000000 0.0000000
#> V2  0.16048825 0.4811226 0.0000000
#> V3 -0.07019999 0.1485488 0.2858419

# Explicit temporal structure
B <- matrix(c(0.4,  0.0, 0.1,
              0.2,  0.3, 0.0,
              0.0, -0.1, 0.5), nrow = 3, byrow = TRUE)
r2 <- simulate_longitudinal(n = 50, tp = 100, vars = 3, temporal = B, seed = 1)

# ESM structure with days and beeps
r3 <- simulate_longitudinal(n = 40, tp = 70, vars = 4,
                             beeps_per_day = 7, seed = 7)
head(r3$data)  # has day + beep columns
#>   id day beep         V1         V2         V3         V4
#> 1  1   1    1  0.5680978  0.9114328 -2.2195456 -2.3066352
#> 2  1   1    2  0.6722628  1.8719812 -1.9654093  0.1321454
#> 3  1   1    3  0.6629147  1.7280159  0.7235085  0.2023174
#> 4  1   1    4 -0.3923350  3.1467353  2.3565133  1.2055313
#> 5  1   1    5 -0.5550650  1.8820806 -0.2582258  1.3479448
#> 6  1   1    6  0.1905947 -0.8165545 -1.0377736 -1.6080693

# With edge-case injection
r4 <- simulate_longitudinal(n = 30, tp = 50, vars = 3,
                             complexity = c("na", "outliers"), seed = 5)
```
