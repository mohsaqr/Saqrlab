# Simulate TNA Network Objects

Simulate multiple TNA network objects (fitted models) for simulation
studies. This function creates TNA model objects by simulating sequence
data and fitting models. For generating raw sequence data with
associated probabilities, use
[`simulate_tna_datasets`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md).
For group TNA models, use
[`simulate_group_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_group_tna_networks.md).

## Usage

``` r
simulate_tna_networks(
  n_networks = 1,
  n_states = 8,
  n_sequences = 100,
  seq_length = 30,
  model_type = "tna",
  use_learning_states = TRUE,
  categories = NULL,
  seed = NULL,
  verbose = TRUE,
  num_sequences = NULL,
  max_seq_length = NULL,
  learning_categories = NULL,
  ...
)

generate_tna_networks(...)
```

## Arguments

- n_networks:

  Integer. Number of networks to generate. Default: 1.

- n_states:

  Integer. Number of states/actions in each network. Default: 5.

- n_sequences:

  Integer. Number of sequences to simulate per network. Default: 100.

- seq_length:

  Integer. Maximum length of each sequence. Default: 30.

- model_type:

  Character. Type of TNA model to fit: "tna", "ftna", "ctna", "atna".
  Default: "tna".

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

- num_sequences:

  Deprecated. Use `n_sequences` instead.

- max_seq_length:

  Deprecated. Use `seq_length` instead.

- learning_categories:

  Deprecated. Use `categories` instead.

- ...:

  Arguments passed to `simulate_tna_networks`.

## Value

A list of length `n_networks` containing fitted TNA model objects. Each
element is a tna model object (class "tna").

## Details

This function differs from
[`simulate_tna_datasets`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
in that it returns only the fitted model objects, not the underlying
sequence data or generating probabilities. Use this when you need TNA
network objects for simulation studies or method comparisons.

**Random Category Selection**: When `categories = NULL` and
`use_learning_states = TRUE`, a random learning category is selected.
Available categories: metacognitive, cognitive, behavioral, social,
motivational, affective, group_regulation.

## See also

[`simulate_group_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_group_tna_networks.md)
for group TNA models,
[`simulate_tna_datasets`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
for generating complete datasets with probabilities,
[`fit_network_model`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md)
for model fitting,
[`get_learning_states`](https://pak.dynasite.org/Saqrlab/reference/get_learning_states.md)
for learning state verbs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 5 TNA networks with random learning category
nets <- simulate_tna_networks(n_networks = 5, seed = 42)

# Generate networks with specific category
meta_nets <- simulate_tna_networks(
  n_networks = 3,
  n_states = 6,
  categories = "metacognitive",
  seed = 123
)

# Generate filtered TNA networks
ftna_nets <- simulate_tna_networks(
  n_networks = 5,
  model_type = "ftna",
  seed = 456
)

# Use the networks
plot(nets[[1]])
tna::centralities(nets[[1]])

# Old parameter names still work
nets <- simulate_tna_networks(
  n_networks = 3,
  num_sequences = 50,
  max_seq_length = 25,
  seed = 42
)
} # }
```
