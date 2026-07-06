# Simulate Social Network Edge List

Generate a simulated edge list for social network analysis. Creates
random connections between nodes with weights and class assignments.

## Usage

``` r
simulate_edge_list(
  n_nodes = 20,
  n_edges = NULL,
  edge_density = 3,
  n_classes = 3,
  directed = TRUE,
  allow_self_loops = FALSE,
  weight_range = c(0.1, 1),
  names = NULL,
  class_probs = NULL,
  seed = NULL
)
```

## Arguments

- n_nodes:

  Integer. Number of nodes (people) in the network. Default: 20.

- n_edges:

  Integer or NULL. Number of edges to generate. If NULL, calculated as
  `n_nodes * edge_density`. Default: NULL.

- edge_density:

  Numeric. Average number of edges per node when `n_edges` is NULL.
  Default: 3.

- n_classes:

  Integer. Number of classes/groups (2-10). Default: 3.

- directed:

  Logical. Whether edges are directed. Default: TRUE.

- allow_self_loops:

  Logical. Whether to allow self-connections. Default: FALSE.

- weight_range:

  Numeric vector of length 2. Range for edge weights. Default: c(0.1,
  1.0).

- names:

  Character vector or NULL. Custom node names. If NULL, uses names from
  `GLOBAL_NAMES`. Default: NULL.

- class_probs:

  Numeric vector or NULL. Probability of each class. Must sum to 1 and
  have length `n_classes`. If NULL, uniform distribution. Default: NULL.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

## Value

A data frame with columns:

- source:

  Character. Name of the source node.

- target:

  Character. Name of the target node.

- weight:

  Numeric. Edge weight in the specified range.

- class:

  Integer. Class assignment (1 to n_classes).

## Details

The function generates a random social network edge list with:

- Nodes named using diverse global names (or custom names)

- Random edges between nodes

- Weights uniformly distributed in the specified range

- Class assignments based on specified probabilities

For undirected networks, each edge appears once (no duplicate A-B, B-A
pairs).

## See also

[`GLOBAL_NAMES`](https://pak.dynasite.org/Saqrlab/reference/GLOBAL_NAMES.md),
[`get_global_names`](https://pak.dynasite.org/Saqrlab/reference/get_global_names.md)

## Examples

``` r
# Basic usage with defaults
edges <- simulate_edge_list(seed = 42)
head(edges)
#>       source    target weight class
#> 1    Anahera  Shoshana 0.7814     1
#> 2    Anahera    Soraya 0.3799     2
#> 3     Bataar    Esther 0.7585     2
#> 4 Cuauhtemoc      Oyun 0.9325     3
#> 5 Cuauhtemoc Oyunbileg 0.4311     2
#> 6     Eloise    Esther 0.1543     1

# Larger network with 5 classes
edges <- simulate_edge_list(
  n_nodes = 50,
  n_classes = 5,
  edge_density = 4,
  seed = 123
)

# Custom names and specific number of edges
edges <- simulate_edge_list(
  n_nodes = 10,
  n_edges = 30,
  names = c("Alice", "Bob", "Carol", "Dave", "Eve",
            "Frank", "Grace", "Hank", "Ivy", "Jack"),
  n_classes = 2,
  seed = 42
)

# Undirected network with unequal class distribution
edges <- simulate_edge_list(
  n_nodes = 30,
  directed = FALSE,
  n_classes = 3,
  class_probs = c(0.5, 0.3, 0.2),
  seed = 42
)

# View class distribution
table(edges$class)
#> 
#>  1  2  3 
#> 52 21 17 
```
