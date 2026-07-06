# Export Simulation Results to CSV

Write a tidy simulation data.frame to CSV. A list of
[`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
objects (or a single one) is first flattened with
[`tidy_simulation_results`](https://pak.dynasite.org/Saqrlab/reference/tidy_simulation_results.md);
a `data.frame` is written as-is.

## Usage

``` r
export_simulation(x, file)
```

## Arguments

- x:

  A `data.frame`, a single
  [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
  object, or a list of `saqr_sim` objects.

- file:

  Character scalar. Path of the CSV file to write.

## Value

The `file` path, invisibly.

## Examples

``` r
s <- run_scenario("power_ttest", seed = 1)
path <- tempfile(fileext = ".csv")
export_simulation(s, path)
```
