# Get a Scenario Recipe

Retrieve the argument-lists for a named scenario. With `case = NULL` the
full list of argument-lists is returned; with an integer `case` the
single argument-list for that case is returned, ready to
[`do.call()`](https://rdrr.io/r/base/do.call.html) against the
scenario's simulator.

## Usage

``` r
get_scenario(scenario, case = NULL)
```

## Arguments

- scenario:

  Character scalar. The scenario name (see
  [`list_scenarios`](https://pak.dynasite.org/Saqrlab/reference/list_scenarios.md)).

- case:

  Integer scalar or `NULL`. When `NULL` (default) the full list of
  argument-lists is returned. When an integer, the single argument-list
  at that index is returned.

## Value

When `case = NULL`, a list of named argument-lists. When `case` is an
integer, a single named argument-list.

## Examples

``` r
get_scenario("power_ttest")
#> [[1]]
#> [[1]]$n_a
#> [1] 30
#> 
#> [[1]]$n_b
#> [1] 30
#> 
#> [[1]]$mean_a
#> [1] 0
#> 
#> [[1]]$mean_b
#> [1] 0.2
#> 
#> 
#> [[2]]
#> [[2]]$n_a
#> [1] 30
#> 
#> [[2]]$n_b
#> [1] 30
#> 
#> [[2]]$mean_a
#> [1] 0
#> 
#> [[2]]$mean_b
#> [1] 0.5
#> 
#> 
#> [[3]]
#> [[3]]$n_a
#> [1] 30
#> 
#> [[3]]$n_b
#> [1] 30
#> 
#> [[3]]$mean_a
#> [1] 0
#> 
#> [[3]]$mean_b
#> [1] 0.8
#> 
#> 
#> [[4]]
#> [[4]]$n_a
#> [1] 100
#> 
#> [[4]]$n_b
#> [1] 100
#> 
#> [[4]]$mean_a
#> [1] 0
#> 
#> [[4]]$mean_b
#> [1] 0.2
#> 
#> 
#> [[5]]
#> [[5]]$n_a
#> [1] 100
#> 
#> [[5]]$n_b
#> [1] 100
#> 
#> [[5]]$mean_a
#> [1] 0
#> 
#> [[5]]$mean_b
#> [1] 0.5
#> 
#> 
#> [[6]]
#> [[6]]$n_a
#> [1] 100
#> 
#> [[6]]$n_b
#> [1] 100
#> 
#> [[6]]$mean_a
#> [1] 0
#> 
#> [[6]]$mean_b
#> [1] 0.8
#> 
#> 
get_scenario("power_ttest", case = 1)
#> $n_a
#> [1] 30
#> 
#> $n_b
#> [1] 30
#> 
#> $mean_a
#> [1] 0
#> 
#> $mean_b
#> [1] 0.2
#> 
```
