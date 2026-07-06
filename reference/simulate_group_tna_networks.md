# Simulate Group TNA Network Objects

Simulate group TNA network objects (fitted group_tna models) for
simulation studies with grouped/clustered data. Each network contains
multiple groups (e.g., classrooms, teams) with their own transition
patterns.

## Usage

``` r
simulate_group_tna_networks(
  n_groups = 5,
  n_actors = 10,
  n_states = 8,
  seq_length_range = c(10, 30),
  use_learning_states = TRUE,
  categories = NULL,
  seed = NULL,
  verbose = TRUE,
  actors_per_group = NULL,
  min_seq_length = NULL,
  max_seq_length = NULL,
  learning_categories = NULL,
  ...
)

generate_group_tna_networks(...)
```

## Arguments

- n_groups:

  Integer. Number of groups in the network. Default: 5.

- n_actors:

  Integer or integer vector. Number of actors per group. Accepts: single
  integer (fixed size), two integers like `c(8, 12)` (min/max), or a
  range like `5:15` (min/max taken from range). Default: 10.

- n_states:

  Integer. Number of states/actions in the network. Default: 5.

- seq_length_range:

  Integer vector of length 2. Range for sequence lengths per actor (min,
  max). Default: c(10, 30).

- use_learning_states:

  Logical. If TRUE, uses learning action verbs as state names. Default:
  TRUE.

- categories:

  Character vector or NULL. Which categories of learning verbs to use.
  Options: "metacognitive", "cognitive", "behavioral", "social",
  "motivational", "affective", "group_regulation", or "all". If NULL
  (default), randomly selects one category.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

- verbose:

  Logical. If TRUE, prints progress messages. Default: TRUE.

- actors_per_group:

  Deprecated. Use `n_actors` instead.

- min_seq_length:

  Deprecated. Use `seq_length_range` instead.

- max_seq_length:

  Deprecated. Use `seq_length_range` instead.

- learning_categories:

  Deprecated. Use `categories` instead.

- ...:

  Arguments passed to `simulate_group_tna_networks`.

## Value

A group_tna model object (class "group_tna") containing:

- networks:

  List of TNA networks, one per group

- data:

  The underlying wide-format data

- group:

  The grouping variable name

## Details

This function generates a single group TNA network with multiple groups.
It simulates long-format sequence data with group structure, converts it
to wide format, and fits a group TNA model using
[`tna::group_model()`](http://sonsoles.me/tna/reference/group_model.md).

**Use Cases**:

- Simulating classroom-level learning behavior data

- Generating team collaboration sequences

- Creating multi-group datasets for method comparison

## See also

[`simulate_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
for individual TNA models,
[`simulate_long_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_long_data.md)
for generating long-format group data,
[`fit_network_model`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md)
for model fitting.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a group TNA network with 4 groups, 15 actors each
group_net <- simulate_group_tna_networks(
  n_groups = 4,
  n_actors = 15,
  n_states = 5,
  seed = 42
)

# Variable group sizes using range notation
var_net <- simulate_group_tna_networks(
  n_groups = 5,
  n_actors = c(8, 15),
  seq_length_range = c(5, 25),
  seed = 123
)

# With specific learning category
ssrl_net <- simulate_group_tna_networks(
  n_groups = 6,
  n_actors = c(10, 20),
  categories = "group_regulation",
  seed = 456
)

# Access individual group networks
names(group_net)
group_net[[1]]$weights

# Old parameter names still work
group_net <- simulate_group_tna_networks(
  actors_per_group = 10,
  min_seq_length = 5,
  max_seq_length = 20,
  seed = 42
)
} # }
```
