# Run a Whole Scenario

Execute every case of a scenario by
[`do.call()`](https://rdrr.io/r/base/do.call.html)-ing its argument-list
against the scenario's simulator. Returns a named list of the resulting
[`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
objects, labelled `case1`, `case2`, ...

## Usage

``` r
run_scenario(scenario, seed = NULL)
```

## Arguments

- scenario:

  Character scalar. The scenario name (see
  [`list_scenarios`](https://pak.dynasite.org/Saqrlab/reference/list_scenarios.md)).

- seed:

  Integer scalar or `NULL`. When supplied, each case is given a
  deterministic seed of `seed + case index` so the whole run is
  reproducible while every case still differs. When `NULL` the
  simulators use their own defaults.

## Value

A named list of
[`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
objects, one per case, named `case1`, `case2`, ...

## Examples

``` r
s <- run_scenario("power_ttest", seed = 1)
s$case1
#> saqr_sim [ttest]  60 x 2  (seed=2)
#>   params: mean_a, mean_b, sd_a, sd_b, n_a, n_b, cohens_d 
#>   cols:   group, score 
```
