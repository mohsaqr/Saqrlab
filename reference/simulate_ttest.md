# Simulate Two-Group Comparison Data (t-test)

Generate data for a two-group comparison with fully specified
ground-truth means and standard deviations. Designed so that
`t.test(score ~ group, data = r$data)` recovers the true difference at
large `n`.

## Usage

``` r
simulate_ttest(
  n_a,
  n_b,
  mean_a,
  mean_b,
  sd_a = 1,
  sd_b = sd_a,
  labels = c("A", "B"),
  seed = NULL
)
```

## Arguments

- n_a:

  Integer. Sample size for group A.

- n_b:

  Integer. Sample size for group B.

- mean_a:

  Numeric. Population mean for group A.

- mean_b:

  Numeric. Population mean for group B.

- sd_a:

  Positive numeric. Standard deviation for group A. Default: 1.

- sd_b:

  Positive numeric. Standard deviation for group B. Default: same as
  `sd_a`.

- labels:

  Character vector of length 2. Group labels. Default: `c("A", "B")`.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `group` (factor) and `score` (numeric).

- `$params`:

  list with `mean_a`, `mean_b`, `sd_a`, `sd_b`, `n_a`, `n_b`, `cohens_d`
  (true Cohen's d based on pooled SD).

## Examples

``` r
r <- simulate_ttest(n_a = 50, n_b = 50, mean_a = 100, mean_b = 105, seed = 1)
t.test(score ~ group, data = r$data)
#> 
#>  Welch Two Sample t-test
#> 
#> data:  score by group
#> t = -27.787, df = 95.793, p-value < 2.2e-16
#> alternative hypothesis: true difference in means between group A and group B is not equal to 0
#> 95 percent confidence interval:
#>  -5.375269 -4.658487
#> sample estimates:
#> mean in group A mean in group B 
#>        100.1004        105.1173 
#> 
r$params$cohens_d
#> [1] 5
```
