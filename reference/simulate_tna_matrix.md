# Simulate TNA Transition Matrix with Node Groupings

Simulate a transition matrix with node groupings, compatible with
`tna::plot_htna()` (hierarchical) and `tna::plot_mlna()` (multilevel)
visualizations. By default creates a 25-node matrix (5 nodes x 5 types)
using learning category names.

This is a convenience wrapper around
[`simulate_htna`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md).

## Usage

``` r
simulate_tna_matrix(
  nodes_per_group = 5,
  group_names = c("Metacognitive", "Cognitive", "Behavioral", "Social", "Motivational"),
  n_groups = 5,
  edge_prob_range = c(0, 1),
  self_loops = FALSE,
  use_learning_states = TRUE,
  categories = c("metacognitive", "cognitive", "behavioral", "social", "motivational"),
  within_prob = 0.4,
  between_prob = 0.15,
  node_prefix = "N",
  seed = NULL,
  verbose = TRUE,
  learning_categories = NULL
)

generate_tna_matrix(...)
```

## Arguments

- nodes_per_group:

  Integer. Number of nodes per group. Default: 5.

- group_names:

  Character vector. Names for each group. Default uses learning
  categories: "Metacognitive", "Cognitive", "Behavioral", "Social",
  "Motivational".

- n_groups:

  Integer. Number of groups. Default: 5.

- edge_prob_range:

  Numeric vector of length 2. Range for edge weights `c(min, max)`.
  Default: c(0, 1).

- self_loops:

  Logical. Allow self-loops (diagonal elements). Default: FALSE.

- use_learning_states:

  Logical. Use learning state verbs as node names. Default: TRUE.

- categories:

  Character vector. Categories for node names, one per group. Default:
  c("metacognitive", "cognitive", "behavioral", "social",
  "motivational").

- within_prob:

  Numeric. Probability of edges within each group. Default: 0.4.

- between_prob:

  Numeric. Probability of edges between groups. Default: 0.15.

- node_prefix:

  Character. Prefix for node names when not using learning states.
  Default: "N".

- seed:

  Integer or NULL. Random seed. Default: NULL.

- verbose:

  Logical. Print progress messages. Default: TRUE.

- learning_categories:

  Deprecated. Use `categories` instead.

- ...:

  Arguments passed to `simulate_tna_matrix`.

## Value

A list with two elements:

- matrix:

  Square transition matrix (rows sum to 1) with named rows/columns.

- node_types:

  Named list mapping group names to node names. Use as `node_types` for
  `plot_htna()` or as `layers` for `plot_mlna()`.

## Details

This function generates a random transition matrix and node groupings.
The output can be used with both hierarchical and multilevel TNA plots:

- For `plot_htna()`: use `net$node_types` directly

- For `plot_mlna()`: use `layers = net$node_types`

## See also

[`simulate_htna`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
for the underlying function,
[`simulate_matrix`](https://pak.dynasite.org/Saqrlab/reference/simulate_matrix.md)
for basic matrix simulation,
[`simulate_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
for TNA model objects.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: 5 groups x 5 nodes = 25 node matrix
net <- simulate_tna_matrix(seed = 42)
net$matrix
net$node_types  # Metacognitive, Cognitive, Behavioral, Social, Motivational

# Use with plot_htna
plot_htna(net$matrix, net$node_types, layout = "polygon")

# Use with plot_mlna
plot_mlna(net$matrix, layers = net$node_types)

# Custom group names (3 groups)
net <- simulate_tna_matrix(
  nodes_per_group = 6,
  group_names = c("Macro", "Meso", "Micro"),
  seed = 42
)

# Custom categories per group
net <- simulate_tna_matrix(
  nodes_per_group = 4,
  group_names = c("Teacher", "Student", "System"),
  categories = c("metacognitive", "cognitive", "behavioral"),
  seed = 123
)
} # }
```
