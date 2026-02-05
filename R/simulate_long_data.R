#' Simulate Long Format Sequence Data
#'
#' @description
#' Generate simulated sequence data in long format, matching the structure of
#' `group_regulation_long` in the tna package. Creates realistic educational
#' data with actors nested in groups, groups nested in courses.
#'
#' @param n_groups Integer. Number of groups. Default: 5.
#' @param n_actors Integer or integer vector of length 2. If single
#'   integer, exact number of actors per group. If vector c(min, max), random
#'   sizes in that range. Default: 10.
#' @param n_courses Integer or character vector. Number of courses or specific
#'   course names (e.g., c("A", "B", "C")). Groups are distributed evenly
#'   across courses. Default: 3.
#' @param n_states Integer. Number of actions to sample when using categories.
#'   Ignored if `states` is provided. Default: 9.
#' @param states Character vector. The action/state names to use. If NULL,
#'   uses learning states based on `categories`. Default: NULL.
#' @param use_learning_states Logical. If TRUE and `states` is NULL, uses
#'   learning state verbs from the specified `categories`. Default: TRUE.
#' @param categories Character vector. Categories of learning states
#'   to use when `states = NULL`. Options: "metacognitive", "cognitive",
#'   "behavioral", "social", "motivational", "affective", "group_regulation",
#'   or "all". Default: "group_regulation".
#' @param seq_length_range Integer vector of length 2. Range for sequence
#'   lengths per actor (min, max). Default: c(10, 30).
#' @param achiever_levels Character vector. Levels for the Achiever variable.
#'   Default: c("High", "Low").
#' @param achiever_probs Numeric vector. Probabilities for each achiever level.
#'   Must sum to 1. Default: c(0.5, 0.5).
#' @param start_time POSIXct or character. Start time for timestamps.
#'   Default: "2025-01-01 10:00:00".
#' @param time_interval_range Numeric vector of length 2. Range for time
#'   intervals between actions in seconds (min, max). Default: c(60, 600).
#' @param trans_matrix Matrix or NULL. Custom transition probability matrix.
#'   If NULL, generates random probabilities. Default: NULL.
#' @param init_probs Numeric vector or NULL. Custom initial probabilities.
#'   If NULL, generates random probabilities. Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @param actors_per_group Deprecated. Use `n_actors` instead.
#' @param actions Deprecated. Use `states` instead.
#' @param action_categories Deprecated. Use `categories` instead.
#' @param n_actions Deprecated. Use `n_states` instead.
#' @param transition_probs Deprecated. Use `trans_matrix` instead.
#' @param initial_probs Deprecated. Use `init_probs` instead.
#'
#' @return A tibble (data frame) with columns:
#' \describe{
#'   \item{Actor}{Integer. Actor identifier (1 to n_actors).}
#'   \item{Achiever}{Character. Achievement level (e.g., "High", "Low").}
#'   \item{Group}{Numeric. Group identifier (1 to n_groups).}
#'   \item{Course}{Character. Course identifier.}
#'   \item{Time}{POSIXct. Timestamp of the action.}
#'   \item{Action}{Character. The action/state performed.}
#' }
#'
#' @details
#' The data structure follows a nested hierarchy:
#' \itemize{
#'   \item **Courses** contain multiple **Groups**
#'   \item **Groups** contain multiple **Actors** (e.g., 10 per group)
#'   \item **Actors** have multiple **Actions** over time
#' }
#'
#' This matches the structure of `group_regulation_long` from the tna package:
#' \itemize{
#'   \item Each group has a fixed number of actors
#'   \item Groups are distributed across courses
#'   \item Actors within a group share the same course
#' }
#'
#' **Built-in Group Regulation Actions**: When `categories = "group_regulation"`,
#' uses the 9 SSRL (Socially Shared Regulation of Learning) actions:
#' adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis.
#'
#' @examples
#' \dontrun{
#' # Basic usage: 5 groups, 10 actors each, 3 courses (new defaults)
#' data <- simulate_long_data(seed = 42)
#'
#' # Explicit new parameter names
#' data <- simulate_long_data(
#'   n_groups = 10,
#'   n_actors = 15,
#'   n_states = 6,
#'   categories = c("metacognitive", "cognitive"),
#'   seq_length_range = c(15, 40),
#'   seed = 42
#' )
#'
#' # Check structure
#' length(unique(data$Actor))   # 50 actors
#' length(unique(data$Group))   # 5 groups
#' table(table(data$Actor, data$Group) > 0)  # 10 actors per group
#'
#' # Variable group sizes
#' data <- simulate_long_data(
#'   n_groups = 30,
#'   n_actors = c(8, 12),
#'   states = c("Read", "Write", "Discuss", "Plan", "Review"),
#'   seed = 456
#' )
#'
#' # Using specific learning state categories
#' data <- simulate_long_data(
#'   n_groups = 50,
#'   n_actors = 10,
#'   categories = c("metacognitive", "cognitive"),
#'   n_states = 8,
#'   seed = 789
#' )
#'
#' # Old parameter names still work (backward compatible)
#' data <- simulate_long_data(
#'   actors_per_group = 10,
#'   actions = c("A", "B", "C"),
#'   seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_onehot_data}} for one-hot encoded format,
#' \code{\link{simulate_sequences}} for generating wide-format sequences,
#' \code{\link{get_learning_states}} for available learning states,
#' \code{\link{generate_tna_networks}} for generating TNA network objects.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @importFrom stats runif
#' @export
simulate_long_data <- function(n_groups = 5,
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
                                # Backward compatibility - old parameter names
                                actors_per_group = NULL,
                                actions = NULL,
                                action_categories = NULL,
                                n_actions = NULL,
                                transition_probs = NULL,
                                initial_probs = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(actors_per_group)) n_actors <- actors_per_group
  if (!is.null(actions)) states <- actions
  if (!is.null(action_categories)) categories <- action_categories
  if (!is.null(n_actions)) n_states <- n_actions
  if (!is.null(transition_probs)) trans_matrix <- transition_probs
  if (!is.null(initial_probs)) init_probs <- initial_probs

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Input Validation ---
  stopifnot(
    is.numeric(n_groups), length(n_groups) == 1, n_groups >= 1,
    is.numeric(n_actors), length(n_actors) %in% c(1, 2),
    all(n_actors >= 1),
    is.numeric(seq_length_range), length(seq_length_range) == 2,
    seq_length_range[1] <= seq_length_range[2], seq_length_range[1] >= 1,
    is.character(achiever_levels), length(achiever_levels) >= 1,
    is.numeric(achiever_probs), length(achiever_probs) == length(achiever_levels),
    abs(sum(achiever_probs) - 1) < 1e-6,
    is.numeric(time_interval_range), length(time_interval_range) == 2,
    time_interval_range[1] <= time_interval_range[2], time_interval_range[1] >= 0
  )

  # Normalize n_actors to range
  if (length(n_actors) == 1) {
    n_actors <- c(n_actors, n_actors)
  }

  # --- Setup Course Names ---
  if (is.numeric(n_courses) && length(n_courses) == 1) {
    course_names <- LETTERS[1:n_courses]
  } else if (is.character(n_courses)) {
    course_names <- n_courses
  } else {
    stop("n_courses must be a single integer or a character vector of course names.")
  }
  n_courses_num <- length(course_names)

  # --- Determine Actions/States ---
  if (is.null(states)) {
    if (use_learning_states) {
      if ("group_regulation" %in% categories) {
        # Built-in group regulation actions (SSRL)
        group_reg_actions <- c(
          "adapt", "cohesion", "consensus", "coregulate",
          "discuss", "emotion", "monitor", "plan", "synthesis"
        )
        if (length(categories) == 1 && categories == "group_regulation") {
          states <- group_reg_actions
        } else {
          # Combine with other categories
          other_cats <- setdiff(categories, "group_regulation")
          other_actions <- if (length(other_cats) > 0) {
            get_learning_states(other_cats, n = max(1, n_states - length(group_reg_actions)))
          } else {
            character(0)
          }
          states <- unique(c(group_reg_actions, other_actions))
        }
      } else {
        # Use learning states from specified categories
        states <- get_learning_states(categories, n = n_states)
      }
    } else {
      # Default to letters
      if (n_states <= 26) {
        states <- LETTERS[1:n_states]
      } else {
        states <- paste0("S", 1:n_states)
      }
    }
  }

  n_action_types <- length(states)
  stopifnot(n_action_types >= 2)

  # --- Generate Transition Probabilities ---
  if (is.null(trans_matrix)) {
    trans_matrix <- seqHMM::simulate_transition_probs(n_action_types)
  }
  dimnames(trans_matrix) <- list(states, states)

  if (is.null(init_probs)) {
    init_probs <- seqHMM::simulate_initial_probs(n_action_types)
  }
  names(init_probs) <- states

  # --- Parse Start Time ---
  if (is.character(start_time)) {
    start_time <- as.POSIXct(start_time)
  }

  # --- Assign Groups to Courses ---
  # Distribute groups across courses evenly
  groups_per_course <- ceiling(n_groups / n_courses_num)
  group_to_course <- rep(course_names, each = groups_per_course)[1:n_groups]

  # --- Create Groups and Actors ---
  all_rows <- list()
  actor_id <- 0

  for (group_id in 1:n_groups) {
    # Determine number of actors in this group
    if (n_actors[1] == n_actors[2]) {
      n_actors_in_group <- n_actors[1]
    } else {
      n_actors_in_group <- sample(n_actors[1]:n_actors[2], 1)
    }

    course <- group_to_course[group_id]

    # Create actors for this group
    for (a in 1:n_actors_in_group) {
      actor_id <- actor_id + 1

      # Assign achiever level
      achiever <- sample(achiever_levels, 1, prob = achiever_probs)

      # Determine sequence length for this actor
      seq_len <- sample(seq_length_range[1]:seq_length_range[2], 1)

      # Generate sequence using Markov chain
      sequence <- character(seq_len)
      sequence[1] <- sample(states, 1, prob = init_probs)

      if (seq_len > 1) {
        for (j in 2:seq_len) {
          current_state <- sequence[j - 1]
          current_idx <- match(current_state, states)
          sequence[j] <- sample(states, 1, prob = trans_matrix[current_idx, ])
        }
      }

      # Generate timestamps
      time_intervals <- runif(seq_len, time_interval_range[1], time_interval_range[2])
      # Add some randomness to start time per actor
      actor_start <- start_time + runif(1, 0, 3600) # Random offset up to 1 hour
      timestamps <- actor_start + cumsum(c(0, time_intervals[-1]))

      # Create data frame for this actor
      actor_df <- data.frame(
        Actor = actor_id,
        Achiever = achiever,
        Group = group_id,
        Course = course,
        Time = timestamps,
        Action = sequence,
        stringsAsFactors = FALSE
      )

      all_rows[[length(all_rows) + 1]] <- actor_df
    }
  }

  # --- Combine All Data ---
  result <- do.call(rbind, all_rows)

  # Sort by Actor, then Time
  result <- result[order(result$Actor, result$Time), ]
  rownames(result) <- NULL

  # Convert to tibble for consistency with tna package
  result <- structure(
    result,
    class = c("tbl_df", "tbl", "data.frame")
  )

  return(result)
}


#' Group Regulation Actions
#'
#' @description
#' The 9 standard group regulation (SSRL) actions used in collaborative
#' learning research.
#'
#' @format A character vector of 9 actions.
#'
#' @details
#' These actions represent Socially Shared Regulation of Learning (SSRL):
#' \describe{
#'   \item{adapt}{Adjusting strategies based on feedback}
#'   \item{cohesion}{Building and maintaining group unity}
#'   \item{consensus}{Reaching group agreement}
#'   \item{coregulate}{Mutually regulating each other's learning}
#'   \item{discuss}{Engaging in group discussion}
#'   \item{emotion}{Expressing or regulating emotions}
#'   \item{monitor}{Tracking and evaluating progress}
#'   \item{plan}{Planning group activities and strategies}
#'   \item{synthesis}{Integrating and combining ideas}
#' }
#'
#' @examples
#' GROUP_REGULATION_ACTIONS
#'
#' # Use in simulation
#' data <- simulate_long_data(
#'   n_groups = 10,
#'   states = GROUP_REGULATION_ACTIONS,
#'   seed = 42
#' )
#'
#' @export
GROUP_REGULATION_ACTIONS <- c(
  "adapt",
  "cohesion",
  "consensus",
  "coregulate",
  "discuss",
  "emotion",
  "monitor",
  "plan",
  "synthesis"
)
