# Simulate Network Matrix

Generate a simple transition matrix for TNA/Markov network analysis.
Each call randomly selects a learning category for node names.

## Usage

``` r
simulate_matrix(
  n_nodes = 9,
  matrix_type = c("transition", "frequency", "co-occurrence", "adjacency"),
  edge_prob = 0.3,
  weighted = TRUE,
  weight_range = c(0, 1),
  directed = TRUE,
  allow_self_loops = FALSE,
  names = NULL,
  seed = NULL
)
```

## Arguments

- n_nodes:

  Integer. Number of nodes. Default: 9.

- matrix_type:

  Character. Type of matrix to generate:

  - `"transition"`: Directed, row-normalized (rows sum to 1)

  - `"frequency"`: Directed, integer counts

  - `"co-occurrence"`: Symmetric

  - `"adjacency"`: Binary or weighted edges

  Default: "transition".

- edge_prob:

  Numeric. Probability of edges existing. Default: 0.3.

- weighted:

  Logical. If TRUE, generates weighted edges. Default: TRUE.

- weight_range:

  Numeric vector of length 2. Range for edge weights. Default: c(0, 1).

- directed:

  Logical. If TRUE, generates directed matrix. Default: TRUE.

- allow_self_loops:

  Logical. If TRUE, allows diagonal entries. Default: FALSE.

- names:

  Character vector or NULL. Custom node names. If NULL, randomly selects
  learning states from a random category. Default: NULL.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

## Value

A numeric matrix with row and column names set to learning states.

## Details

This function generates a simple network matrix. Each call randomly
picks one learning category (metacognitive, cognitive, behavioral,
social, motivational, affective, or group_regulation) and uses verbs
from that category as node names.

For matrices with multiple node types, use
[`simulate_htna`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md).

## See also

[`simulate_htna`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
for multi-type matrices

## Examples

``` r
# Simple 9-node transition matrix
mat <- simulate_matrix(seed = 42)
mat
#>            Regulate   Plan  Judge Reflect Monitor Forecast Anticipate  Check
#> Regulate     0.0000 0.2400 0.0000  0.0000  0.2466   0.0000          0 0.5134
#> Plan         1.0000 0.0000 0.0000  0.0000  0.0000   0.0000          0 0.0000
#> Judge        0.0000 0.4135 0.0000  0.1031  0.0000   0.0000          0 0.0000
#> Reflect      0.0000 0.0000 0.0000  0.0000  0.0000   0.0000          0 0.0000
#> Monitor      0.0000 0.1702 0.3583  0.3466  0.0000   0.1249          0 0.0000
#> Forecast     0.5343 0.0000 0.0000  0.2381  0.0000   0.0000          0 0.0377
#> Anticipate   0.0000 0.3534 0.1845  0.0000  0.2857   0.1764          0 0.0000
#> Check        0.0000 0.1316 0.4461  0.0000  0.3049   0.1174          0 0.0000
#> Adapt        0.0000 0.0000 0.0000  0.0000  0.7792   0.2208          0 0.0000
#>             Adapt
#> Regulate   0.0000
#> Plan       0.0000
#> Judge      0.4833
#> Reflect    0.0000
#> Monitor    0.0000
#> Forecast   0.1899
#> Anticipate 0.0000
#> Check      0.0000
#> Adapt      0.0000
rowSums(mat)  # Rows sum to 1
#>   Regulate       Plan      Judge    Reflect    Monitor   Forecast Anticipate 
#>     1.0000     1.0000     0.9999     0.0000     1.0000     1.0000     1.0000 
#>      Check      Adapt 
#>     1.0000     1.0000 

# Frequency matrix
mat <- simulate_matrix(n_nodes = 5, matrix_type = "frequency", seed = 42)

# Co-occurrence matrix (symmetric)
mat <- simulate_matrix(n_nodes = 6, matrix_type = "co-occurrence", seed = 42)
isSymmetric(mat)  # TRUE
#> [1] TRUE

# Custom names
mat <- simulate_matrix(n_nodes = 4, names = c("A", "B", "C", "D"), seed = 42)
```
