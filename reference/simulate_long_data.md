# Simulate Long Format Sequence Data

Generate simulated sequence data in long format, matching the structure
of `group_regulation_long` in the tna package. Creates realistic
educational data with actors nested in groups, groups nested in courses.

## Usage

``` r
simulate_long_data(
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

- Action:

  Character. The action/state performed.

## Details

The data structure follows a nested hierarchy:

- **Courses** contain multiple **Groups**

- **Groups** contain multiple **Actors** (e.g., 10 per group)

- **Actors** have multiple **Actions** over time

This matches the structure of `group_regulation_long` from the tna
package:

- Each group has a fixed number of actors

- Groups are distributed across courses

- Actors within a group share the same course

**Built-in Group Regulation Actions**: When
`categories = "group_regulation"`, uses the 9 SSRL (Socially Shared
Regulation of Learning) actions: adapt, cohesion, consensus, coregulate,
discuss, emotion, monitor, plan, synthesis.

## See also

[`simulate_onehot_data`](https://pak.dynasite.org/Saqrlab/reference/simulate_onehot_data.md)
for one-hot encoded format,
[`simulate_sequences`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
for generating wide-format sequences,
[`get_learning_states`](https://pak.dynasite.org/Saqrlab/reference/get_learning_states.md)
for available learning states,
[`simulate_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
for generating TNA network objects.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage: 5 groups, 10 actors each, 3 courses (new defaults)
data <- simulate_long_data(seed = 42)

# Explicit new parameter names
data <- simulate_long_data(
  n_groups = 10,
  n_actors = 15,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seq_length_range = c(15, 40),
  seed = 42
)

# Check structure
length(unique(data$Actor))   # 50 actors
length(unique(data$Group))   # 5 groups
table(table(data$Actor, data$Group) > 0)  # 10 actors per group

# Variable group sizes
data <- simulate_long_data(
  n_groups = 30,
  n_actors = c(8, 12),
  states = c("Read", "Write", "Discuss", "Plan", "Review"),
  seed = 456
)

# Using specific learning state categories
data <- simulate_long_data(
  n_groups = 50,
  n_actors = 10,
  categories = c("metacognitive", "cognitive"),
  n_states = 8,
  seed = 789
)

# Old parameter names still work (backward compatible)
data <- simulate_long_data(
  actors_per_group = 10,
  actions = c("A", "B", "C"),
  seed = 42
)
} # }
```
