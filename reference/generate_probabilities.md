# Generate Transition and Initial Probabilities

Generate random transition probabilities and initial state probabilities
for a Markov chain using the seqHMM package.

## Usage

``` r
generate_probabilities(
  n_states = 8,
  states = NULL,
  alpha = 1,
  diag_c = 0,
  seed = NULL,
  possible_state_names = NULL
)
```

## Arguments

- n_states:

  Integer. Number of states in the Markov chain. Default: 5.

- states:

  Character vector. Names for the states. Must have at least `n_states`
  elements. If NULL, uses letters A, B, C, ... Default: NULL.

- alpha:

  Numeric. Dirichlet concentration parameter passed to
  [`seqHMM::simulate_initial_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html)
  and
  [`seqHMM::simulate_transition_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html).
  Small values (e.g., 0.1) produce sparse matrices with a few dominant
  transitions; large values (e.g., 10) produce near-uniform matrices.
  Default: 1.

- diag_c:

  Numeric. Diagonal boost added before row normalisation, passed to
  [`seqHMM::simulate_transition_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html).
  Higher values create "sticky" states with strong self-transitions.
  Default: 0.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

- possible_state_names:

  Deprecated. Use `states` instead.

## Value

A list containing:

- initial_probs:

  Named numeric vector of initial state probabilities.

- transition_probs:

  Square matrix of transition probabilities with row and column names
  set to state names.

- state_names:

  Character vector of state names used.

## Details

The transition probabilities are generated using
[`seqHMM::simulate_transition_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html)
and initial probabilities using
[`seqHMM::simulate_initial_probs()`](https://rdrr.io/pkg/seqHMM/man/simulate_pars.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate probabilities for 4 states with default names
probs <- generate_probabilities(n_states = 4, seed = 42)

# Generate probabilities with custom state names
probs <- generate_probabilities(
  n_states = 4,
  states = c("A", "B", "C", "D", "E"),
  seed = 123
)

# View initial probabilities
probs$initial_probs

# View transition matrix
probs$transition_probs

# Old parameter name still works
probs <- generate_probabilities(
  n_states = 3,
  possible_state_names = c("X", "Y", "Z")
)
} # }
```
