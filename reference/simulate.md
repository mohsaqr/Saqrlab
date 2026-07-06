# Simulate Data via a Unified Dispatcher

A single discoverable entry point over every explicit-parameter
simulator in the package. You pass the simulation `type` plus the named
arguments that simulator expects, and `simulate()` forwards them
straight through, returning that simulator's
[`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object unchanged. Use
[`list_simulators`](https://pak.dynasite.org/Saqrlab/reference/list_simulators.md)
to see the available types.

For random, seed-driven generation (where the seed picks both the data
and the structural parameters) and a bare `data.frame` return, use the
separate
[`simulate_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
entry point instead.

## Usage

``` r
simulate(type, ..., seed = NULL)
```

## Arguments

- type:

  Character scalar. The simulation type. One of the values in the `type`
  column of
  [`list_simulators`](https://pak.dynasite.org/Saqrlab/reference/list_simulators.md):
  `"ttest"`, `"anova"`, `"correlation"`, `"clusters"`, `"prediction"`,
  `"regression"`, `"lpa"`, `"lca"`, `"fa"`, `"seq_clusters"`,
  `"longitudinal"`, `"mlm"`, `"growth"`, `"irt"`, `"survival"`, `"hmm"`.

- ...:

  Named arguments passed straight through to the matching simulator. See
  the help page of the underlying simulator (e.g.
  [`simulate_ttest`](https://pak.dynasite.org/Saqrlab/reference/simulate_ttest.md))
  for the arguments it accepts.

- seed:

  Integer or NULL. Random seed, forwarded to the target simulator.

## Value

The [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object returned by the dispatched simulator (with `$data`, `$params`,
`$type`, `$seed`).

## See also

[`list_simulators`](https://pak.dynasite.org/Saqrlab/reference/list_simulators.md)
for the catalogue of types,
[`simulate_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
for random seed-driven generation.

## Examples

``` r
simulate("ttest", n_a = 50, n_b = 50, mean_a = 0, mean_b = 0.5, seed = 1)
#> saqr_sim [ttest]  100 x 2  (seed=1)
#>   params: mean_a, mean_b, sd_a, sd_b, n_a, n_b, cohens_d 
#>   cols:   group, score 

simulate("anova", n = 30, means = c(10, 12, 15), seed = 1)
#> saqr_sim [anova]  90 x 2  (seed=1)
#>   params: means, sds, n, labels, eta_squared 
#>   cols:   group, score 
```
