# Score Parameter Recovery Against Simulated Ground Truth

Compare a method's point estimates to the known true values that
generated a simulation, matched **by name**. This is the core "did my
method recover the truth I simulated?" check of the simulation
laboratory: feed in a
[`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object (or a plain named list / vector of true values) and a named
vector of estimates, and get back a tidy one-row-per-parameter table of
errors and tolerance flags.

## Usage

``` r
validate_recovery(
  sim,
  estimates,
  params = NULL,
  tolerance = 0.1,
  relative = TRUE
)

# S3 method for class 'recovery_result'
print(x, ...)

# S3 method for class 'recovery_result'
summary(object, ...)
```

## Arguments

- sim:

  A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
  object (its `$params` supplies the ground truth) *or* a plain named
  list / named numeric vector of true values.

- estimates:

  A named numeric vector (or named list of scalars) of point estimates
  produced by a fitted method. Names are matched against the (flattened)
  names of the truth.

- params:

  Optional character vector restricting which parameter names to
  compare. Default `NULL` uses the intersection of names available in
  both the truth and the estimates.

- tolerance:

  Positive numeric. Threshold a parameter must fall within to count as
  recovered. Interpreted as a relative tolerance when `relative = TRUE`,
  otherwise an absolute tolerance. Default `0.1`.

- relative:

  Logical. If `TRUE` (default), `within_tol` compares
  `rel_error <= tolerance`; if `FALSE`, compares
  `abs_error <= tolerance`.

- x:

  A `recovery_result` object.

- ...:

  Further arguments (ignored).

- object:

  A `recovery_result` object.

## Value

An S3 object of class `c("recovery_result", "data.frame")` with one row
per compared parameter and columns:

- `parameter`:

  Parameter name (flattened, e.g. `coefs.x1`).

- `true`:

  True value from the simulation.

- `estimate`:

  Estimated value.

- `abs_error`:

  `abs(estimate - true)`.

- `rel_error`:

  `abs_error / abs(true)`, `NA` when `true == 0`.

- `within_tol`:

  Logical recovery flag (see `relative`).

Use `summary.recovery_result` for a one-row scorecard.

## Examples

``` r
# t-test: recover the true Cohen's d
sim <- simulate_ttest(n_a = 200, n_b = 200, mean_a = 0, mean_b = 0.5,
                      sd_a = 1, sd_b = 1, seed = 1)
fit <- t.test(score ~ group, data = sim$data)
d_hat <- unname(fit$estimate[2] - fit$estimate[1])   # mean_b - mean_a
validate_recovery(sim, estimates = c(cohens_d = d_hat))
#> Parameter recovery  (relative tolerance = 0.1)
#>   1 parameters | 1 within tolerance (100.0%)
#>   mean abs error = 0.005098 | mean rel error = 0.0102
#> 
#>  parameter true  estimate   abs_error rel_error within_tol
#>   cohens_d  0.5 0.5050981 0.005098052 0.0101961       TRUE

# Regression: recover known coefficients from lm()
rsim <- simulate_regression(
  coefs = c("(Intercept)" = 1, x1 = 2, x2 = -1.5),
  predictor_sds = c(x1 = 1, x2 = 1), error_sd = 1, n = 2000, seed = 7
)
fit_lm <- lm(y ~ x1 + x2, data = rsim$data)
est <- stats::setNames(coef(fit_lm),
                       paste0("coefs.", names(coef(fit_lm))))
validate_recovery(rsim, estimates = est)
#> Parameter recovery  (relative tolerance = 0.1)
#>   3 parameters | 3 within tolerance (100.0%)
#>   mean abs error = 0.01819 | mean rel error = 0.012
#> 
#>          parameter true  estimate   abs_error   rel_error within_tol
#>  coefs.(Intercept)  1.0  0.994209 0.005790962 0.005790962       TRUE
#>           coefs.x1  2.0  1.986102 0.013897903 0.006948951       TRUE
#>           coefs.x2 -1.5 -1.465106 0.034894381 0.023262921       TRUE
```
