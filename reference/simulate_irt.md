# Simulate Item Response Theory (IRT) Data

Generate item-response data under a known IRT model with ground-truth
item and person parameters. Supports the 1PL (Rasch), 2PL, 3PL
dichotomous models and the graded response model (GRM) for ordered
polytomous items. Designed so that, at large `n_persons`, item
proportion-correct recovers item difficulty (negatively) and person
total score recovers true ability (positively).

Response probabilities are computed for the whole person-by-item grid in
one vectorised pass ([`outer()`](https://rdrr.io/r/base/outer.html) /
matrix ops) and responses are drawn with a single comparison against a
uniform random matrix – no loops.

## Usage

``` r
simulate_irt(
  n_persons = 500,
  n_items = 20,
  model = c("2PL", "1PL", "3PL", "GRM"),
  n_categories = 2,
  a = NULL,
  b = NULL,
  c = NULL,
  theta = NULL,
  seed = NULL
)
```

## Arguments

- n_persons:

  Integer. Number of persons (rows). Default: 500.

- n_items:

  Integer. Number of items (columns). Default: 20.

- model:

  Character. One of `"2PL"`, `"1PL"`, `"3PL"`, `"GRM"`. Matched via
  [`match.arg`](https://rdrr.io/r/base/match.arg.html).

- n_categories:

  Integer. Number of response categories. Must be 2 for the dichotomous
  models (1PL/2PL/3PL) and `> 2` for GRM. Default: 2.

- a:

  Numeric vector of length `n_items` or NULL. Item discriminations. NULL
  auto-generates `rlnorm(meanlog = 0, sdlog = 0.3)` (values near 1). For
  1PL, discrimination is fixed at 1 regardless of `a`.

- b:

  Numeric vector of length `n_items` or NULL. Item difficulties. NULL
  auto-generates `N(0, 1)`. For GRM, `b` is the item location around
  which ordered category thresholds are centred.

- c:

  Numeric vector of length `n_items` or NULL. Lower-asymptote (guessing)
  parameters. NULL auto-generates `0.2` for 3PL and `0` otherwise.
  Ignored (forced to 0) for 1PL/2PL/GRM.

- theta:

  Numeric vector of length `n_persons` or NULL. True person abilities.
  NULL auto-generates `N(0, 1)`.

- seed:

  Integer or NULL. Random seed.
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called only when
  not NULL.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with `n_persons` rows and columns `item1`...`itemJ`. Values
  are `0/1` for dichotomous models and `0`...`n_categories - 1` for GRM.

- `$params`:

  list with `model`, `a`, `b`, `c` (as used), `theta` (true abilities),
  `thresholds` (GRM category thresholds matrix, NULL otherwise), and
  `n_categories`.

## Examples

``` r
# 2PL with default auto-generated parameters
fit <- simulate_irt(n_persons = 500, n_items = 20, model = "2PL", seed = 1)
print(fit)
#> saqr_sim [irt]  500 x 20  (seed=1)
#>   params: model, a, b, c, theta, thresholds, n_categories 
#>   cols:   item1, item2, item3, item4, item5, item6, ..., item20 
fit$params$b
#>       item1       item2       item3       item4       item5       item6 
#> -0.73732753  0.29066665 -0.88484957  0.20800648 -0.04773017 -1.68452065 
#>       item7       item8       item9      item10      item11      item12 
#> -0.14422656  1.18021367  0.68139992  0.14324763 -1.19231644  1.16922865 
#>      item13      item14      item15      item16      item17      item18 
#>  0.07920171 -0.45177375  1.64202821 -0.76959232  0.30336096  1.28173742 
#>      item19      item20 
#>  0.60222280 -0.30702226 

# Recover difficulty: easier items are answered correctly more often
prop_correct <- colMeans(as.data.frame(fit))
cor(prop_correct, fit$params$b)        # strongly negative
#> [1] -0.9664313

# Recover ability: higher total score tracks higher true theta
total <- rowSums(as.data.frame(fit))
cor(total, fit$params$theta)           # strongly positive
#> [1] 0.870793

# Rasch (1PL): discrimination fixed at 1
simulate_irt(n_persons = 300, n_items = 10, model = "1PL", seed = 2)
#> saqr_sim [irt]  300 x 10  (seed=2)
#>   params: model, a, b, c, theta, thresholds, n_categories 
#>   cols:   item1, item2, item3, item4, item5, item6, ..., item10 

# Graded response model with 4 ordered categories
simulate_irt(n_persons = 400, n_items = 8, model = "GRM",
             n_categories = 4, seed = 3)
#> saqr_sim [irt]  400 x 8  (seed=3)
#>   params: model, a, b, c, theta, thresholds, n_categories 
#>   cols:   item1, item2, item3, item4, item5, item6, item7, item8 
```
