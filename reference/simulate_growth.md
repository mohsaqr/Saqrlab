# Simulate Latent Growth-Curve Data

Generate longitudinal data from a linear latent growth-curve model. Each
subject draws a random intercept and slope from a bivariate normal with
specified means, standard deviations, and correlation, then \$\$y\_{it}
= b\_{0i} + b\_{1i} \\ t + e\_{it},\$\$ with `time` coded
`0..(n_time - 1)`. The result is returned in long format. Designed so
that the mean per-subject slope recovers `slope_mean` and the spread of
per-subject intercepts recovers `intercept_sd` at large `n`.

## Usage

``` r
simulate_growth(
  n = 200,
  n_time = 5,
  intercept_mean = 0,
  slope_mean = 1,
  intercept_sd = 1,
  slope_sd = 0.5,
  intercept_slope_cor = 0,
  residual_sd = 1,
  seed = NULL
)
```

## Arguments

- n:

  Integer. Number of subjects. Default: 200.

- n_time:

  Integer \>= 2. Number of measurement occasions per subject. Default:
  5.

- intercept_mean:

  Numeric. Mean of subject intercepts \\b\_{0i}\\. Default: 0.

- slope_mean:

  Numeric. Mean of subject slopes \\b\_{1i}\\. Default: 1.

- intercept_sd:

  Positive numeric. Standard deviation of subject intercepts. Default:
  1.

- slope_sd:

  Positive numeric. Standard deviation of subject slopes. Default: 0.5.

- intercept_slope_cor:

  Numeric in `[-1, 1]`. Correlation between subject intercepts and
  slopes. Default: 0.

- residual_sd:

  Positive numeric. Level-1 (occasion) residual standard deviation
  \\\sigma_e\\. Default: 1.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  long-format data.frame with columns `subject` (factor), `time`
  (numeric, `0..n_time-1`), and `y` (numeric).

- `$params`:

  list with `means` (intercept/slope), `sds` (intercept/slope),
  `correlation`, `residual_sd`, `n_time`, `subject_intercepts` (true
  \\b\_{0i}\\), and `subject_slopes` (true \\b\_{1i}\\).

## Examples

``` r
r <- simulate_growth(n = 300, n_time = 6, slope_mean = 2,
                     intercept_sd = 1.5, seed = 1)
print(r)
#> saqr_sim [growth]  1800 x 3  (seed=1)
#>   params: means, sds, correlation, residual_sd, n_time, subject_intercepts, subject_slopes 
#>   cols:   subject, time, y 
head(r)
#>   subject time         y
#> 1       1    0 -1.280748
#> 2       1    1  3.009581
#> 3       1    2  4.482301
#> 4       1    3  6.943021
#> 5       1    4  8.710993
#> 6       1    5 10.157770

# Correlated intercepts and slopes
r2 <- simulate_growth(n = 250, intercept_slope_cor = 0.5, seed = 42)
r2$params$correlation
#> [1] 0.5
```
