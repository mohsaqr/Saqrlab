# Simulate Prediction/Regression Data with Known Coefficients

Generate a regression dataset with both continuous and categorical
predictors, a non-linear term, and known ground-truth coefficients. A
richer version of
[`simulate_regression`](https://pak.dynasite.org/Saqrlab/reference/simulate_regression.md)
for testing prediction workflows.

## Usage

``` r
simulate_prediction(
  n,
  coefs,
  cat_levels = NULL,
  cat_effects = NULL,
  error_sd = 1,
  predictor_means = NULL,
  predictor_sds = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Integer. Sample size.

- coefs:

  Named numeric vector. Must include names matching continuous
  predictors (e.g. `x1`, `x2`). May include `"(Intercept)"`. Non-linear
  terms are NOT automatically generated; this specifies the linear part.

- cat_levels:

  Named list. Each element is a character vector of factor levels for a
  categorical predictor. Default: `NULL` (no categorical predictors).

- cat_effects:

  Named list of numeric vectors. Effect of each level (same length and
  names as `cat_levels`). Default: random effects.

- error_sd:

  Positive numeric. Residual standard deviation. Default: 1.

- predictor_means:

  Named numeric vector or NULL. Means for continuous predictors.
  Default: all zeros.

- predictor_sds:

  Named numeric vector or NULL. SDs for continuous predictors. Default:
  all ones.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `y`, continuous predictors, and categorical
  predictors (factors).

- `$params`:

  list with `coefs`, `cat_effects`, `error_sd`, `predictor_means`,
  `predictor_sds`, `r_squared` (population \\R^2\\).

## Examples

``` r
r <- simulate_prediction(
  n = 200,
  coefs = c("(Intercept)" = 5, x1 = 2, x2 = -1),
  cat_levels = list(treatment = c("control", "drug_a", "drug_b")),
  cat_effects = list(treatment = c(0, 3, 5)),
  error_sd = 2, seed = 42
)
summary(lm(y ~ ., data = r$data))
#> 
#> Call:
#> lm(formula = y ~ ., data = r$data)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -5.6238 -1.2751  0.0521  1.3570  4.5064 
#> 
#> Coefficients:
#>                 Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)       4.9896     0.2391  20.866  < 2e-16 ***
#> x1                1.8857     0.1389  13.578  < 2e-16 ***
#> x2               -1.0033     0.1447  -6.933 5.89e-11 ***
#> treatmentdrug_a   3.0066     0.3356   8.958 2.65e-16 ***
#> treatmentdrug_b   4.6081     0.3322  13.870  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.897 on 195 degrees of freedom
#> Multiple R-squared:  0.6975, Adjusted R-squared:  0.6913 
#> F-statistic: 112.4 on 4 and 195 DF,  p-value: < 2.2e-16
#> 
r$params$r_squared
#> [1] 0.6933315
```
