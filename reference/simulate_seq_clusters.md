# Simulate Sequence Data with Known Cluster Structure

Generate wide-format sequence data where each row is drawn from one of K
Markov chains (clusters), each governed by its own transition matrix.
Supports explicit matrix supply or automatic random generation.

## Usage

``` r
simulate_seq_clusters(
  trans_list = NULL,
  props = NULL,
  n = 300L,
  seq_length = 20L,
  init_probs = NULL,
  n_clusters = 3L,
  n_states = 10L,
  states = NULL,
  seed = NULL
)
```

## Arguments

- trans_list:

  Named or unnamed list of square numeric transition matrices, or `NULL`
  for automatic generation. When `NULL`, `n_clusters` and `n_states` are
  used to build random row-stochastic matrices.

- props:

  Numeric vector of mixing proportions (need not sum to 1; normalised
  internally). Length must equal `K`. Defaults to equal mixing.

- n:

  Positive integer. Total number of sequences to generate.

- seq_length:

  Positive integer. Number of time-points per sequence (columns
  `T1`...`T{seq_length}`). Default 20.

- init_probs:

  Either a numeric vector of initial state probabilities shared across
  clusters, a list of per-cluster vectors, or `NULL` (uniform). Must be
  named consistently with the state names derived from the transition
  matrix row/column names.

- n_clusters:

  Positive integer. Number of clusters when `trans_list = NULL`. Default
  3.

- n_states:

  Positive integer. Number of states when `trans_list = NULL`. Default
  10.

- states:

  Character vector of state names when `trans_list = NULL`. Defaults to
  `paste0("S", seq_len(n_states))`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

## Value

A named list with elements:

- `data`:

  data.frame with columns `T1`...`T{seq_length}` (character state
  labels) and `true_cluster` (integer 1...K).

- `params`:

  list with `trans_list` (the K matrices used), `props` (normalised
  mixing proportions), and `init_probs` (per-cluster initial probability
  vectors as a list).

## Examples

``` r
m1 <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE,
             dimnames = list(c("A","B"), c("A","B")))
m2 <- matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2, byrow = TRUE,
             dimnames = list(c("A","B"), c("A","B")))
r  <- simulate_seq_clusters(trans_list = list(m1, m2), n = 100, seed = 1)
head(r$data)
#>   T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20
#> 1  A  B  A  A  B  A  B  A  A   B   B   B   A   B   A   B   B   A   B   A
#> 2  A  B  A  B  B  A  B  A  B   A   B   A   B   A   A   B   A   B   B   A
#> 3  A  A  A  A  A  A  A  A  A   B   B   B   B   B   B   B   B   B   B   B
#> 4  B  A  A  A  B  B  B  B  A   A   A   B   A   A   A   B   B   A   A   B
#> 5  B  A  A  B  B  A  B  B  B   A   B   A   B   B   A   B   A   A   B   B
#> 6  B  B  B  B  B  B  B  B  B   A   B   B   A   B   A   A   A   B   A   A
#>   true_cluster
#> 1            2
#> 2            2
#> 3            1
#> 4            1
#> 5            2
#> 6            1
r$params$props
#> [1] 0.5 0.5
```
