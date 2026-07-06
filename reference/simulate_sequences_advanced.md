# Simulate Markov Chain Sequences (Advanced)

Generate sequences of states from a Markov chain with advanced stability
and instability modes. Supports stable transitions, probability
perturbation, and unlikely jump modes.

## Usage

``` r
simulate_sequences_advanced(
  n_sequences = 1000,
  seq_length = 20,
  n_states = 8,
  states = NULL,
  use_learning_states = TRUE,
  categories = "all",
  trans_matrix = NULL,
  init_probs = NULL,
  stable_transitions = NULL,
  stability_prob = 0.95,
  unstable_mode = "unlikely_jump",
  unstable_random_transition_prob = 0.4,
  unstable_perturb_noise = 0.5,
  unlikely_prob_threshold = 0.1,
  na_range = c(0, 0),
  include_na = TRUE,
  seed = NULL,
  transition_matrix = NULL,
  initial_probabilities = NULL,
  num_rows = NULL,
  max_seq_length = NULL,
  min_na = NULL,
  max_na = NULL
)
```

## Arguments

- n_sequences:

  Integer. Number of sequences (rows) to generate. Default: 1000.

- seq_length:

  Integer. Maximum length of each sequence. Default: 20.

- n_states:

  Integer. Number of states when auto-generating probabilities. Ignored
  if `trans_matrix` is provided. Default: 5.

- states:

  Character vector. Names for states when auto-generating. If NULL, uses
  letters (A, B, C, ...) or learning states if enabled. Ignored if
  `trans_matrix` is provided. Default: NULL.

- use_learning_states:

  Logical. If TRUE and auto-generating, uses realistic learning action
  verbs as state names. Default: TRUE.

- categories:

  Character vector. Categories of learning states to use. Only used if
  `use_learning_states = TRUE`. Default: "all".

- trans_matrix:

  Square numeric matrix of transition probabilities. Rows must sum to 1.
  Row names define state names. If NULL, random probabilities are
  generated. Default: NULL.

- init_probs:

  Named numeric vector of initial state probabilities. Must sum to 1. If
  NULL, random probabilities are generated. Default: NULL.

- stable_transitions:

  List of character vectors. Each vector contains two state names
  defining a stable transition pair (from, to). Default: NULL (no stable
  transitions).

- stability_prob:

  Numeric in (0 to 1). Probability of following a stable transition when
  in a stable state. Default: 0.95.

- unstable_mode:

  Character. Mode for unstable transitions. One of:

  - "random_jump": Uniform random jump to any state.

  - "perturb_prob": Perturb transition probabilities with noise.

  - "unlikely_jump": Jump to states with low probability.

  Default: "unlikely_jump".

- unstable_random_transition_prob:

  Numeric in (0 to 1). Probability of taking an unstable action.
  Default: 0.4.

- unstable_perturb_noise:

  Numeric in (0 to 1). Noise factor for probability perturbation mode.
  Default: 0.5.

- unlikely_prob_threshold:

  Numeric in (0 to 1). Threshold below which transitions are considered
  "unlikely". Default: 0.1.

- na_range:

  Integer vector of length 2 (min, max) or single integer (min=max).
  Range of NA values per sequence. Default: c(0, 0).

- include_na:

  Logical. Whether to include NAs in sequences. Default: TRUE.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

- transition_matrix:

  Deprecated. Use `trans_matrix` instead.

- initial_probabilities:

  Deprecated. Use `init_probs` instead.

- num_rows:

  Deprecated. Use `n_sequences` instead.

- max_seq_length:

  Deprecated. Use `seq_length` instead.

- min_na:

  Deprecated. Use `na_range` instead.

- max_na:

  Deprecated. Use `na_range` instead.

## Value

A data frame with `n_sequences` rows and `seq_length` columns. Each row
is a sequence of state names, potentially with trailing NAs.

## Details

The function extends basic Markov chain simulation with:

**Stable Transitions**: If the current state has a defined stable
transition and a random draw is below `stability_prob`, the sequence
follows the stable transition directly.

**Unstable Modes** (when not following stable transitions):

- "random_jump": With probability `unstable_random_transition_prob`,
  jump uniformly to any state.

- "perturb_prob": Multiply transition probabilities by random noise in
  range `[1-noise, 1+noise]`, then renormalize.

- "unlikely_jump": With probability `unstable_random_transition_prob`,
  jump to a state with transition probability below threshold.

NAs are added at the end of sequences, preserving at least 2 non-NA
values.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simplest usage: all defaults
sequences <- simulate_sequences_advanced(seed = 42)

# Create a 4-state transition matrix
trans_mat <- matrix(c(
  0.6, 0.2, 0.1, 0.1,
  0.1, 0.7, 0.1, 0.1,
  0.1, 0.1, 0.6, 0.2,
  0.2, 0.1, 0.1, 0.6
), nrow = 4, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C", "D")

init_probs <- c(A = 0.25, B = 0.25, C = 0.25, D = 0.25)

# Define stable transitions: A->B and C->D are "stable"
stable <- list(c("A", "B"), c("C", "D"))

# Generate sequences with unlikely_jump mode
sequences <- simulate_sequences_advanced(
  trans_matrix = trans_mat,
  init_probs = init_probs,
  seq_length = 30,
  n_sequences = 100,
  stable_transitions = stable,
  stability_prob = 0.95,
  unstable_mode = "unlikely_jump",
  unstable_random_transition_prob = 0.3,
  na_range = c(0, 5)
)

# Old parameter names still work (backward compatible)
sequences <- simulate_sequences_advanced(
  transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 30,
  num_rows = 100
)
} # }
```
