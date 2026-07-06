# Simulate TNA Datasets (Sequences + Models + Probabilities)

Simulate complete TNA datasets including simulated sequences, fitted
models, and generating probabilities. This function combines the full
simulation workflow: probability generation -\> sequence simulation -\>
model fitting.

Use
[`simulate_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
if you only need the fitted models.

Supports using realistic learning action verbs as state names for
educational research simulations.

## Usage

``` r
simulate_tna_datasets(
  n_datasets = 1,
  n_states = 8,
  n_sequences = 100,
  seq_length = 30,
  states = NULL,
  use_learning_states = TRUE,
  categories = NULL,
  smart_select = TRUE,
  na_range = c(0, 0),
  model_type = "tna",
  use_advanced = FALSE,
  stable_transitions = NULL,
  stability_prob = 0.95,
  unstable_mode = "random_jump",
  unstable_random_transition_prob = 0.5,
  include_data = TRUE,
  include_probs = TRUE,
  seed = NULL,
  verbose = TRUE,
  state_names = NULL,
  learning_categories = NULL,
  num_rows = NULL,
  max_seq_length = NULL,
  min_na = NULL,
  max_na = NULL,
  include_probabilities = NULL
)

generate_tna_datasets(...)

generate_sequence_data(...)
```

## Arguments

- n_datasets:

  Integer. Number of datasets to generate. Default: 1.

- n_states:

  Integer. Number of states in each network. Default: 5.

- n_sequences:

  Integer. Number of sequences to simulate per network. Default: 100.

- seq_length:

  Integer. Maximum length of each sequence. Default: 30.

- states:

  Character vector. Names for the states. If NULL and
  `use_learning_states = FALSE`, uses letters (A, B, C, ...). Ignored if
  `use_learning_states = TRUE`. Default: NULL.

- use_learning_states:

  Logical. If TRUE, uses learning action verbs as state names (e.g.,
  "Plan", "Monitor", "Read"). Default: TRUE.

- categories:

  Character vector. Which categories of learning verbs to use. Options:
  "metacognitive", "cognitive", "behavioral", "social", "motivational",
  "affective", "group_regulation", or "all". If NULL and
  `use_learning_states = TRUE`, a random category is selected. Only used
  if `use_learning_states = TRUE`. Default: NULL.

- smart_select:

  Logical. If TRUE, intelligently selects learning states based on
  n_states (small networks use fewer categories). Only used if
  `use_learning_states = TRUE`. Default: TRUE.

- na_range:

  Integer vector of length 2 (min, max) or single integer. Range of NA
  values per sequence. Default: c(0, 0).

- model_type:

  Character. Type of TNA model to fit. One of: "tna", "ftna", "ctna",
  "atna". Default: "tna".

- use_advanced:

  Logical. If TRUE, uses
  [`simulate_sequences_advanced()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences_advanced.md)
  with stability modes. If FALSE, uses basic
  [`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md).
  Default: FALSE.

- stable_transitions:

  List of character vectors defining stable transitions (only used if
  `use_advanced = TRUE`). Default: NULL.

- stability_prob:

  Numeric (0 to 1). Probability of following stable transitions (only
  used if `use_advanced = TRUE`). Default: 0.95.

- unstable_mode:

  Character. Mode for unstable transitions: "random_jump",
  "perturb_prob", or "unlikely_jump". (Only used if
  `use_advanced = TRUE`). Default: "random_jump".

- unstable_random_transition_prob:

  Numeric (0 to 1). Probability of unstable action (only used if
  `use_advanced = TRUE`). Default: 0.5.

- include_data:

  Logical. If TRUE, includes the generating sequence data in the output.
  Default: TRUE.

- include_probs:

  Logical. If TRUE, includes the generating transition matrix and
  initial probabilities in the output. Default: TRUE.

- seed:

  Integer or NULL. Random seed for reproducibility. If NULL, no seed is
  set. Default: NULL.

- verbose:

  Logical. If TRUE, prints progress messages. Default: TRUE.

- state_names:

  Deprecated. Use `states` instead.

- learning_categories:

  Deprecated. Use `categories` instead.

- num_rows:

  Deprecated. Use `n_sequences` instead.

- max_seq_length:

  Deprecated. Use `seq_length` instead.

- min_na:

  Deprecated. Use `na_range` instead.

- max_na:

  Deprecated. Use `na_range` instead.

- include_probabilities:

  Deprecated. Use `include_probs` instead.

- ...:

  Arguments passed to `simulate_tna_datasets`.

## Value

A list of length `n_datasets`. Each element is a list containing:

- model:

  The fitted TNA model object.

- transition_probs:

  The generating transition matrix (if `include_probs = TRUE`).

- initial_probs:

  The generating initial probabilities (if `include_probs = TRUE`).

- sequences:

  The simulated sequence data frame (if `include_data = TRUE`).

- params:

  List of parameters used for this network.

## Details

The function generates sequence data by:

1.  Generating random transition and initial probabilities using
    [`seqHMM::simulate_transition_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html)
    and
    [`seqHMM::simulate_initial_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html).

2.  Simulating sequences from those probabilities using either
    [`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
    or
    [`simulate_sequences_advanced()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences_advanced.md).

3.  Fitting a TNA model to the sequences using
    [`fit_network_model()`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md).

**Learning States**: When `use_learning_states = TRUE`, state names are
drawn from a curated collection of 180+ student learning action verbs
organized into 7 categories:

- **metacognitive**: Plan, Monitor, Evaluate, Reflect, Regulate, ...

- **cognitive**: Read, Study, Analyze, Summarize, Memorize, ...

- **behavioral**: Practice, Annotate, Research, Review, Write, ...

- **social**: Collaborate, Discuss, Explain, Share, Teach, ...

- **motivational**: Focus, Persist, Explore, Strive, Commit, ...

- **affective**: Enjoy, Appreciate, Cope, Manage, Curious, ...

- **group_regulation**: Adapt, Cohesion, Consensus, Coregulate, ...

## See also

[`simulate_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
for generating TNA network objects,
[`get_learning_states`](https://pak.dynasite.org/Saqrlab/reference/get_learning_states.md)
for retrieving learning state verbs,
[`select_states`](https://pak.dynasite.org/Saqrlab/reference/select_states.md)
for intelligent state selection,
[`list_learning_categories`](https://pak.dynasite.org/Saqrlab/reference/list_learning_categories.md)
for viewing available categories,
[`simulate_sequences`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
for basic sequence simulation,
[`fit_network_model`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md)
for model fitting.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simplest usage - generates 1 dataset with sequences + model + probabilities
data <- simulate_tna_datasets(seed = 42)
data[[1]]$sequences        # The sequence data
data[[1]]$model            # The fitted TNA model
data[[1]]$transition_probs # The generating probabilities

# Generate 5 datasets with learning states
data <- simulate_tna_datasets(
  n_datasets = 5,
  n_states = 4,
  seed = 42
)

# Generate with specific learning category
learning_data <- simulate_tna_datasets(
  n_datasets = 5,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

# View the state names
learning_data[[1]]$params$state_names
# e.g., c("Plan", "Monitor", "Read", "Practice", "Discuss", "Focus")

# Generate with letter names (disable learning states)
letter_data <- simulate_tna_datasets(
  n_datasets = 5,
  n_states = 8,
  use_learning_states = FALSE,
  seed = 123
)
} # }
```
