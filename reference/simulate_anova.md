# Simulate Multi-Group Comparison Data (ANOVA)

Generate data for a one-way ANOVA with fully specified group means and
standard deviations. Designed so that
`aov(score ~ group, data = r$data)` recovers the true group effects at
large `n`.

## Usage

``` r
simulate_anova(n, means, sds = 1, labels = NULL, seed = NULL)
```

## Arguments

- n:

  Integer or integer vector. If a single value, the per-group sample
  size (equal groups). If a vector, must have length equal to
  `length(means)` giving per-group sizes.

- means:

  Numeric vector. Population mean for each group.

- sds:

  Numeric scalar or vector. Standard deviation(s) for each group. Scalar
  is recycled. Default: 1.

- labels:

  Character vector or NULL. Group labels. Default:
  `paste0("G", seq_along(means))`.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `group` (factor) and `score` (numeric).

- `$params`:

  list with `means`, `sds` (as vector), `n` (per-group sizes), `labels`,
  `eta_squared` (true population \\\eta^2\\).

## Examples

``` r
r <- simulate_anova(n = 30, means = c(10, 12, 15), seed = 1)
summary(aov(score ~ group, data = r$data))
#>             Df Sum Sq Mean Sq F value Pr(>F)    
#> group        2  383.5   191.7   238.8 <2e-16 ***
#> Residuals   87   69.9     0.8                   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
r$params$eta_squared
#> [1] 0.8137045

# Unequal groups and heteroscedastic
r2 <- simulate_anova(n = c(50, 30, 20), means = c(5, 5, 10),
                      sds = c(1, 2, 3), seed = 42)
```
