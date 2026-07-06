# List Available Simulators

Return the catalogue of simulation types reachable through
[`simulate`](https://pak.dynasite.org/Saqrlab/reference/simulate.md),
one row per type, with the underlying function name and a one-line
description. Print it directly.

Note: this lists the explicit-parameter simulators (those returning a
[`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)).
The separate
[`simulate_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
entry point provides random, seed-driven generation and returns a bare
`data.frame`; it is not part of this catalogue.

## Usage

``` r
list_simulators()
```

## Value

A base `data.frame` with one row per simulator and columns:

- `type`:

  The `type` string passed to
  [`simulate`](https://pak.dynasite.org/Saqrlab/reference/simulate.md).

- `function`:

  The name of the underlying simulator function.

- `description`:

  A short one-line description.

## See also

[`simulate`](https://pak.dynasite.org/Saqrlab/reference/simulate.md)

## Examples

``` r
list_simulators()
#>            type              function
#> 1         ttest        simulate_ttest
#> 2         anova        simulate_anova
#> 3   correlation  simulate_correlation
#> 4      clusters     simulate_clusters
#> 5    prediction   simulate_prediction
#> 6    regression   simulate_regression
#> 7           lpa          simulate_lpa
#> 8           lca          simulate_lca
#> 9            fa           simulate_fa
#> 10 seq_clusters simulate_seq_clusters
#> 11 longitudinal simulate_longitudinal
#> 12          mlm          simulate_mlm
#> 13       growth       simulate_growth
#> 14          irt          simulate_irt
#> 15     survival     simulate_survival
#> 16          hmm          simulate_hmm
#>                                                                         description
#> 1                               Two-group comparison with known means/SDs (t-test).
#> 2                    Multi-group comparison with known group means (one-way ANOVA).
#> 3          Multivariate normal data from an explicit correlation/covariance matrix.
#> 4                              Gaussian mixture with known cluster centres and SDs.
#> 5  Regression data with continuous + categorical predictors and known coefficients.
#> 6                   Linear regression data with known coefficients and residual SD.
#> 7               Latent profile analysis: continuous indicators from known profiles.
#> 8          Latent class analysis: binary indicators from known class probabilities.
#> 9             Factor analysis: indicators from known loadings and factor structure.
#> 10              Categorical sequences drawn from known cluster-specific generators.
#> 11               Longitudinal panel data with a known within-person VAR(1) process.
#> 12             Multilevel (mixed-effects) data with known fixed and random effects.
#> 13         Latent growth-curve data with known intercept/slope means and variances.
#> 14           Item response data from a known IRT model (difficulty/discrimination).
#> 15            Survival/time-to-event data with known covariate effects (Cox-style).
#> 16                 Hidden Markov sequences from known transition/emission matrices.
```
