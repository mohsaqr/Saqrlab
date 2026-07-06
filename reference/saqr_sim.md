# Create a saqr_sim Object

Wraps simulation output in a standardised S3 class. Every simulation
function that returns numerical data uses this class, giving a
consistent interface: `$data` (the data.frame), `$params` (ground-truth
generating parameters), `$type`, `$seed`.

Backward-compatible: `$data` and `$params` still work exactly as before.
Additionally, `[` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) delegate
to the underlying data, so existing code that subscripts or coerces the
result keeps working.

## Usage

``` r
saqr_sim(data, params, type, seed = NULL, call = NULL, extras = list())
```

## Arguments

- data:

  A data.frame (or matrix) of simulated data.

- params:

  A list of ground-truth generating parameters.

- type:

  Character label identifying the simulation type (e.g. `"lpa"`,
  `"regression"`, `"longitudinal"`).

- seed:

  The random seed actually used (integer or NULL).

- call:

  The matched call that created this object (optional).

- extras:

  Named list of additional fields to attach (optional).

## Value

An object of class `saqr_sim` (inherits from `list`).

## Examples

``` r
sim <- saqr_sim(
  data = data.frame(x = rnorm(10), y = rnorm(10)),
  params = list(mean_x = 0, mean_y = 0),
  type = "demo", seed = 1
)
sim
#> saqr_sim [demo]  10 x 2  (seed=1)
#>   params: mean_x, mean_y 
#>   cols:   x, y 
names(sim)
#> [1] "data"   "params" "type"   "seed"   "call"  
head(sim)
#>             x          y
#> 1  2.02334405  0.0155075
#> 2  0.86249250 -1.6209591
#> 3 -0.02490949 -0.6654647
#> 4  0.60063495 -0.5748405
#> 5  1.21648074 -0.9018930
#> 6 -1.17653155  1.4915994
dim(sim)
#> [1] 10  2
sim$params
#> $mean_x
#> [1] 0
#> 
#> $mean_y
#> [1] 0
#> 
```
