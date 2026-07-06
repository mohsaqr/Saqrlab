# Generate Parameter Grid for Simulations

Generate a grid of parameter combinations for simulation studies using
various sampling methods. Can be called with no arguments for a demo
grid.

## Usage

``` r
generate_param_grid(param_ranges = NULL, n = 10, method = "random")

create_param_grid(...)
```

## Arguments

- param_ranges:

  Named list of parameter ranges. Each element can be:

  - A numeric vector of length 2 (min, max) for continuous/integer
    parameters.

  - A vector of values for categorical parameters.

  If NULL (default), uses demo ranges for TNA simulation:
  `list(n_sequences = c(50, 500), seq_length = c(10, 50), n_states = c(4, 12))`.

- n:

  Integer. Number of parameter combinations to generate. Default: 10.

- method:

  Character. Sampling method. One of:

  "random"

  :   Random uniform sampling within ranges.

  "grid"

  :   Regular grid sampling (may exceed n, then subsampled).

  "lhs"

  :   Latin Hypercube Sampling for better coverage (requires lhs
      package).

  Default: "random".

- ...:

  Arguments passed to `generate_param_grid`.

## Value

A data frame with `n` rows and columns for each parameter.

## Details

**Method Details:**

- "random": Draws uniform random values within each range. Integer
  parameters (detected when min and max are both integers) are rounded.

- "grid": Creates a regular grid with approximately `n^(1/d)` points per
  dimension (where d is the number of parameters). If the resulting grid
  exceeds n points, it is randomly subsampled.

- "lhs": Uses Latin Hypercube Sampling for space-filling designs that
  provide better coverage of the parameter space than random sampling.

Categorical parameters are sampled uniformly with replacement for all
methods.

## Examples

``` r
# Simplest usage: demo grid with default TNA parameters
grid <- generate_param_grid()
head(grid)
#>   n_sequences seq_length n_states
#> 1         490         15        8
#> 2         343         43       10
#> 3         413         21        4
#> 4         460         48       10
#> 5         271         36        6
#> 6          71         16       11

# Define custom parameter ranges
ranges <- list(
  num_rows = c(50, 500),       # Integer parameter
  max_seq_length = c(10, 100), # Integer parameter
  stability_prob = c(0.7, 1.0) # Continuous parameter
)

# Random sampling
grid_random <- generate_param_grid(ranges, n = 20, method = "random")

# Latin Hypercube Sampling
grid_lhs <- generate_param_grid(ranges, n = 20, method = "lhs")

# Grid sampling
grid_regular <- generate_param_grid(ranges, n = 20, method = "grid")
```
