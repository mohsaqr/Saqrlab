# Simulate Cluster Data with Known Centers

Generate multivariate data from a mixture of Gaussians with fully
specified cluster centers and standard deviations. Designed for testing
clustering algorithms where the ground truth is known.

## Usage

``` r
simulate_clusters(n, centers, sds = 1, props = NULL, seed = NULL)
```

## Arguments

- n:

  Integer or integer vector. If a single value, the total sample size
  (allocated by `props`). If a vector of length `k`, the per-cluster
  sizes.

- centers:

  Numeric matrix (`k x d`). Each row is a cluster centroid. `k` = number
  of clusters, `d` = number of dimensions.

- sds:

  Numeric scalar, vector of length `k`, or matrix (`k x d`). Standard
  deviations per cluster (and optionally per dimension). Scalar is
  recycled. Default: 1.

- props:

  Numeric vector of length `k` or NULL. Mixing proportions when `n` is a
  single value. Normalised internally. Default: equal mixing.

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `x1`...`xd` and integer column `true_cluster`
  (1...k).

- `$params`:

  list with `centers` (matrix), `sds` (k x d matrix), `props`
  (normalised), `n` (per-cluster sizes).

## Examples

``` r
centers <- matrix(c(0, 0,
                    5, 5,
                    10, 0), nrow = 3, byrow = TRUE)
r <- simulate_clusters(n = 300, centers = centers, seed = 1)
plot(x2 ~ x1, data = r$data, col = r$data$true_cluster)

r$params$centers
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    5    5
#> [3,]   10    0

# Per-cluster sizes and SDs
r2 <- simulate_clusters(n = c(100, 50, 150), centers = centers,
                         sds = c(0.5, 1.0, 2.0), seed = 42)
```
