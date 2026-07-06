# Simulate Two-Level Multilevel (Hierarchical) Data

Generate data for a two-level model where level-1 units are nested
within clusters. The outcome follows \$\$y\_{ij} = (\beta_0 + u\_{0j}) +
\sum_k \beta_k x\_{kij} + (u\_{1j} x\_{1ij}) + e\_{ij},\$\$ where the
first predictor carries the fixed `slope`, cluster random intercepts
\\u\_{0j} \sim N(0, \tau\_{00})\\ with \\\tau\_{00} = \frac{icc}{1 -
icc}\\\sigma_e^2\\, and (optionally) a cluster random slope \\u\_{1j}
\sim N(0, \mathrm{slope\\sd}^2)\\ on the first predictor. Designed so
that the intraclass correlation and fixed slope are recovered at large
`n_clusters`.

## Usage

``` r
simulate_mlm(
  n_clusters = 30,
  cluster_size = 20,
  intercept = 0,
  slope = 0.5,
  residual_sd = 1,
  icc = 0.1,
  n_predictors = 1,
  random_slope = FALSE,
  slope_sd = 0,
  seed = NULL
)
```

## Arguments

- n_clusters:

  Integer. Number of level-2 clusters. Default: 30.

- cluster_size:

  Integer scalar or integer vector of length `n_clusters`. Number of
  level-1 units per cluster. A vector gives an unbalanced design.
  Default: 20.

- intercept:

  Numeric. Fixed intercept \\\beta_0\\. Default: 0.

- slope:

  Numeric. Fixed slope \\\beta_1\\ on the first predictor (`x1`).
  Default: 0.5.

- residual_sd:

  Positive numeric. Level-1 residual standard deviation \\\sigma_e\\.
  Default: 1.

- icc:

  Numeric in `[0, 1)`. Target intraclass correlation; used to derive the
  random-intercept variance \\\tau\_{00} = \frac{icc}{1 -
  icc}\\\sigma_e^2\\. Default: 0.1.

- n_predictors:

  Integer \>= 1. Number of level-1 predictors (`x1`..`xK`). The first
  carries `slope`; the remaining predictors get fixed coefficients drawn
  from `N(0, 1)`. Default: 1.

- random_slope:

  Logical. If `TRUE`, add a cluster-level random slope on the first
  predictor. Default: `FALSE`.

- slope_sd:

  Non-negative numeric. Standard deviation of the random slope
  \\u\_{1j}\\ (used only when `random_slope = TRUE`). Default: 0.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `cluster_id` (factor), `y` (numeric), and
  `x1`..`xK` (numeric predictors).

- `$params`:

  list with `intercept`, `slope`, `betas` (full fixed-effect vector),
  `icc`, `tau00`, `residual_sd`, `slope_sd`, `random_slope`,
  `cluster_size`, `cluster_intercepts` (true \\u\_{0j}\\), and
  `cluster_slopes` (true \\u\_{1j}\\, `NULL` unless `random_slope`).

## Examples

``` r
r <- simulate_mlm(n_clusters = 40, cluster_size = 25, slope = 0.8,
                  icc = 0.2, seed = 1)
print(r)
#> saqr_sim [mlm]  1000 x 3  (seed=1)
#>   params: intercept, slope, betas, icc, tau00, residual_sd, slope_sd, random_slope, cluster_size, cluster_intercepts, cluster_slopes 
#>   cols:   cluster_id, y, x1 
r$params$tau00
#> [1] 0.25

# Unbalanced clusters with a random slope
r2 <- simulate_mlm(n_clusters = 30, cluster_size = sample(10:30, 30, TRUE),
                   random_slope = TRUE, slope_sd = 0.3, seed = 42)
summary(r2)
#> Simulation type: mlm
#> Seed: 42
#> 
#> --- Data ---
#>    cluster_id        y                   x1          
#>  19     : 30   Min.   :-3.816233   Min.   :-3.01793  
#>  30     : 29   1st Qu.:-0.864847   1st Qu.:-0.66052  
#>  2      : 28   Median :-0.001917   Median :-0.03349  
#>  15     : 28   Mean   :-0.063216   Mean   :-0.02464  
#>  22     : 28   3rd Qu.: 0.638855   3rd Qu.: 0.65320  
#>  1      : 26   Max.   : 3.932679   Max.   : 3.22907  
#>  (Other):432                                         
#> 
#> --- Parameters ---
#> List of 11
#>  $ intercept         : num 0
#>  $ slope             : num 0.5
#>  $ betas             : Named num 0.5
#>   ..- attr(*, "names")= chr "x1"
#>  $ icc               : num 0.1
#>  $ tau00             : num 0.111
#>  $ residual_sd       : num 1
#>  $ slope_sd          : num 0.3
#>  $ random_slope      : logi TRUE
#>  $ cluster_size      : Named int [1:30] 26 28 10 13 19 21 17 24 23 24 ...
#>   ..- attr(*, "names")= chr [1:30] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
#>  $ cluster_intercepts: Named num [1:30] 0.457 -0.188 0.121 0.211 0.135 ...
#>   ..- attr(*, "names")= chr [1:30] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
#>  $ cluster_slopes    : Named num [1:30] 0.137 0.211 0.311 -0.183 0.151 ...
#>   ..- attr(*, "names")= chr [1:30] "cluster_1" "cluster_2" "cluster_3" "cluster_4" ...
```
