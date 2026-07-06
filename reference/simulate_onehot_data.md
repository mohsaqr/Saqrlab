# Simulate One-Hot Encoded Sequence Data

Generate simulated sequence data in one-hot encoded format. Creates the
same hierarchical structure as
[`simulate_long_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_long_data.md)
(Actor, Achiever, Group, Course, Time) but with binary 0/1 columns for
each state instead of a categorical Action column.

## Usage

``` r
simulate_onehot_data(
  n_groups = 5,
  n_actors = 10,
  n_courses = 3,
  n_states = 9,
  states = NULL,
  use_learning_states = TRUE,
  categories = "group_regulation",
  seq_length_range = c(10, 30),
  achiever_levels = c("High", "Low"),
  achiever_probs = c(0.5, 0.5),
  start_time = "2025-01-01 10:00:00",
  time_interval_range = c(60, 600),
  trans_matrix = NULL,
  init_probs = NULL,
  sort_states = FALSE,
  state_prefix = "",
  seed = NULL,
  actors_per_group = NULL,
  actions = NULL,
  action_categories = NULL,
  n_actions = NULL,
  transition_probs = NULL,
  initial_probs = NULL
)
```

## Arguments

- n_groups:

  Integer. Number of groups. Default: 5.

- n_actors:

  Integer or integer vector of length 2. If single integer, exact number
  of actors per group. If vector c(min, max), random sizes in that
  range. Default: 10.

- n_courses:

  Integer or character vector. Number of courses or specific course
  names (e.g., c("A", "B", "C")). Groups are distributed evenly across
  courses. Default: 3.

- n_states:

  Integer. Number of actions to sample when using categories. Ignored if
  `states` is provided. Default: 9.

- states:

  Character vector. The action/state names to use. If NULL, uses
  learning states based on `categories`. Default: NULL.

- use_learning_states:

  Logical. If TRUE and `states` is NULL, uses learning state verbs from
  the specified `categories`. Default: TRUE.

- categories:

  Character vector. Categories of learning states to use when
  `states = NULL`. Options: "metacognitive", "cognitive", "behavioral",
  "social", "motivational", "affective", "group_regulation", or "all".
  Default: "group_regulation".

- seq_length_range:

  Integer vector of length 2. Range for sequence lengths per actor (min,
  max). Default: c(10, 30).

- achiever_levels:

  Character vector. Levels for the Achiever variable. Default: c("High",
  "Low").

- achiever_probs:

  Numeric vector. Probabilities for each achiever level. Must sum to 1.
  Default: c(0.5, 0.5).

- start_time:

  POSIXct or character. Start time for timestamps. Default: "2025-01-01
  10:00:00".

- time_interval_range:

  Numeric vector of length 2. Range for time intervals between actions
  in seconds (min, max). Default: c(60, 600).

- trans_matrix:

  Matrix or NULL. Custom transition probability matrix. If NULL,
  generates random probabilities. Default: NULL.

- init_probs:

  Numeric vector or NULL. Custom initial probabilities. If NULL,
  generates random probabilities. Default: NULL.

- sort_states:

  Logical. Sort state columns alphabetically. Default: FALSE.

- state_prefix:

  Character. Prefix for state column names. Default: "".

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

- actors_per_group:

  Deprecated. Use `n_actors` instead.

- actions:

  Deprecated. Use `states` instead.

- action_categories:

  Deprecated. Use `categories` instead.

- n_actions:

  Deprecated. Use `n_states` instead.

- transition_probs:

  Deprecated. Use `trans_matrix` instead.

- initial_probs:

  Deprecated. Use `init_probs` instead.

## Value

A tibble (data frame) with columns:

- Actor:

  Integer. Actor identifier (1 to n_actors).

- Achiever:

  Character. Achievement level (e.g., "High", "Low").

- Group:

  Numeric. Group identifier (1 to n_groups).

- Course:

  Character. Course identifier.

- Time:

  POSIXct. Timestamp of the action.

- :

  Integer (0/1). One column per state, with 1 indicating the active
  state for that row.

## Details

This function internally calls
[`simulate_long_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_long_data.md)
to generate the base data, then applies an internal one-hot helper to
convert the Action column to one-hot encoded columns.

Each row will have exactly one state column with value 1, and all others
with value 0, representing the action performed at that time point.

## See also

[`simulate_long_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_long_data.md)
for categorical Action format, `Nestimate::action_to_onehot()` for
converting existing long data,
[`simulate_sequences`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
for generating wide-format sequences.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with new standardized params
data <- simulate_onehot_data(n_groups = 5, n_actors = 3, seed = 42)
head(data)

# Custom states with prefix
data <- simulate_onehot_data(
  n_groups = 3,
  states = c("Read", "Write", "Think"),
  state_prefix = "state_",
  seed = 123
)
names(data)

# Verify one-hot encoding (each row sums to 1)
state_cols <- setdiff(names(data), c("Actor", "Achiever", "Group", "Course", "Time"))
all(rowSums(data[, state_cols]) == 1)  # TRUE

# Old parameter names still work (backward compatible)
data <- simulate_onehot_data(
  actors_per_group = 5,
  actions = c("A", "B", "C"),
  seed = 42
)
} # }
```
