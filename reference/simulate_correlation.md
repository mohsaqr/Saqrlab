# Simulate Correlated Multivariate Data

Generate multivariate normal data from an explicit correlation (or
covariance) matrix. Designed so that `cor(r$data)` recovers `sigma` at
large `n`.

## Usage

``` r
simulate_correlation(n, sigma, means = NULL, var_names = NULL, seed = NULL)
```

## Arguments

- n:

  Integer. Sample size.

- sigma:

  Numeric matrix. Either a correlation matrix (all diagonal = 1) or a
  covariance matrix. Must be symmetric and positive-definite (or
  near-PD; corrected internally).

- means:

  Numeric vector or NULL. Population means for each variable. Default:
  all zeros.

- var_names:

  Character vector or NULL. Variable names. Default:
  `paste0("x", seq_len(ncol(sigma)))`.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `x1`...`xp`.

- `$params`:

  list with `sigma` (the input matrix), `means`, `is_correlation`
  (logical).

## Examples

``` r
# Correlation matrix
R <- matrix(c(1, 0.6, 0.3,
              0.6, 1, 0.5,
              0.3, 0.5, 1), nrow = 3)
r <- simulate_correlation(n = 200, sigma = R, seed = 1)
cor(r$data)           # should approximate R
#>           x1        x2        x3
#> x1 1.0000000 0.5556851 0.3162869
#> x2 0.5556851 1.0000000 0.4764833
#> x3 0.3162869 0.4764833 1.0000000
r$params$sigma        # the true matrix
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.6  0.3
#> [2,]  0.6  1.0  0.5
#> [3,]  0.3  0.5  1.0

# With means and custom names
r2 <- simulate_correlation(n = 500, sigma = R,
                            means = c(10, 20, 30),
                            var_names = c("IQ", "GPA", "income"), seed = 42)
```
