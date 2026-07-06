# Simulate Linear Regression Data with Known Coefficients

Generate a dataset suitable for linear regression with fully specified
ground-truth coefficients and predictor standard deviations. Designed so
that `lm(y ~ ., data = r$data)` recovers the true coefficients at large
`n`.

## Usage

``` r
simulate_regression(coefs, predictor_sds, error_sd, n, seed = NULL)
```

## Arguments

- coefs:

  Named numeric vector. May include `"(Intercept)"`. All other names
  become predictor column names. Predictors are generated as independent
  \\N(0, \sigma_j)\\ variables.

- predictor_sds:

  Named numeric vector. Names must exactly match the non-intercept names
  in `coefs`. All values must be positive.

- error_sd:

  Positive numeric. Standard deviation of the residuals.

- n:

  Positive integer. Sample size.

- seed:

  Integer or NULL. Random seed.

## Value

A named list with elements:

- `data`:

  data.frame with column `y` and one column per predictor.

- `params`:

  list with `coefs`, `predictor_sds`, and `error_sd`.

## Examples

``` r
coefs <- c("(Intercept)" = 2, x1 = 3, x2 = -1)
r <- simulate_regression(coefs = coefs,
                          predictor_sds = c(x1 = 1, x2 = 1),
                          error_sd = 0.5, n = 500, seed = 42)
coef(lm(y ~ ., data = r$data))  # should be close to c(2, 3, -1)
#> (Intercept)          x1          x2 
#>   1.9798791   3.0200245  -0.9995026 
```
