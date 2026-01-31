#' Simulate Long Format Sequence Data
#'
#' @description
#' Generate simulated sequence data in long format, matching the structure of
#' `group_regulation_long` in the tna package. Creates realistic educational
#' data with actors nested in groups, groups nested in courses.
#'
#' @param n_groups Integer. Number of groups. Default: 20.
#' @param actors_per_group Integer or integer vector of length 2. If single
#'   integer, exact number of actors per group. If vector c(min, max), random
#'   sizes in that range. Default: 10.
#' @param n_courses Integer or character vector. Number of courses or specific
#'   course names (e.g., c("A", "B", "C")). Groups are distributed evenly
#'   across courses. Default: 3.
#' @param actions Character vector. The action/state names to use. If NULL,
#'   uses learning states based on `action_categories`. Default: NULL.
#' @param action_categories Character vector. Categories of learning states
#'   to use when `actions = NULL`. Options: "metacognitive", "cognitive",
#'   "behavioral", "social", "motivational", "affective", "group_regulation",
#'   or "all". Default: "group_regulation".
#' @param n_actions Integer. Number of actions to sample when using categories.
#'   Ignored if `actions` is provided. Default: 9.
#' @param seq_length_range Integer vector of length 2. Range for sequence
#'   lengths per actor (min, max). Default: c(8, 25).
#' @param achiever_levels Character vector. Levels for the Achiever variable.
#'   Default: c("High", "Low").
#' @param achiever_probs Numeric vector. Probabilities for each achiever level.
#'   Must sum to 1. Default: c(0.5, 0.5).
#' @param start_time POSIXct or character. Start time for timestamps.
#'   Default: "2025-01-01 10:00:00".
#' @param time_interval_range Numeric vector of length 2. Range for time
#'   intervals between actions in seconds (min, max). Default: c(60, 600).
#' @param transition_probs Matrix or NULL. Custom transition probability matrix.
#'   If NULL, generates random probabilities. Default: NULL.
#' @param initial_probs Numeric vector or NULL. Custom initial probabilities.
#'   If NULL, generates random probabilities. Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
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
#' **Built-in Group Regulation Actions**: When `action_categories = "group_regulation"`,
#' uses the 9 SSRL (Socially Shared Regulation of Learning) actions:
#' adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis.
#'
#' @examples
#' \dontrun{
#' # Basic usage: 20 groups, 10 actors each, 3 courses
#' data <- simulate_long_data(
#'   n_groups = 20,
#'   actors_per_group = 10,
#'   seed = 42
#' )
#'
#' # Check structure
#' length(unique(data$Actor))   # 200 actors
#' length(unique(data$Group))   # 20 groups
#' table(table(data$Actor, data$Group) > 0)  # 10 actors per group
#'
#' # Larger dataset like group_regulation_long
#' data <- simulate_long_data(
#'   n_groups = 200,
#'   actors_per_group = 10,
#'   n_courses = 3,
#'   seed = 123
#' )
#'
#' # Custom actions
#' data <- simulate_long_data(
#'   n_groups = 30,
#'   actors_per_group = c(8, 12),
#'   actions = c("Read", "Write", "Discuss", "Plan", "Review"),
#'   seed = 456
#' )
#'
#' # Using learning state categories
#' data <- simulate_long_data(
#'   n_groups = 50,
#'   actors_per_group = 10,
#'   action_categories = c("metacognitive", "cognitive"),
#'   n_actions = 8,
#'   seed = 789
#' )
#'
#' # Check groups per course
#' aggregate(Group ~ Course, data = unique(data[, c("Group", "Course")]),
#'           FUN = length)
#' }
#'
#' @seealso
#' \code{\link{simulate_sequences}} for generating wide-format sequences,
#' \code{\link{get_learning_states}} for available learning states,
#' \code{\link{generate_tna_networks}} for generating TNA network objects.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @importFrom stats runif
#' @export
simulate_long_data <- function(n_groups = 20,
                                actors_per_group = 10,
                                n_courses = 3,
                                actions = NULL,
                                action_categories = "group_regulation",
                                n_actions = 9,
                                seq_length_range = c(8, 25),
                                achiever_levels = c("High", "Low"),
                                achiever_probs = c(0.5, 0.5),
                                start_time = "2025-01-01 10:00:00",
                                time_interval_range = c(60, 600),
                                transition_probs = NULL,
                                initial_probs = NULL,
                                seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Input Validation ---
  stopifnot(
    is.numeric(n_groups), length(n_groups) == 1, n_groups >= 1,
    is.numeric(actors_per_group), length(actors_per_group) %in% c(1, 2),
    all(actors_per_group >= 1),
    is.numeric(seq_length_range), length(seq_length_range) == 2,
    seq_length_range[1] <= seq_length_range[2], seq_length_range[1] >= 1,
    is.character(achiever_levels), length(achiever_levels) >= 1,
    is.numeric(achiever_probs), length(achiever_probs) == length(achiever_levels),
    abs(sum(achiever_probs) - 1) < 1e-6,
    is.numeric(time_interval_range), length(time_interval_range) == 2,
    time_interval_range[1] <= time_interval_range[2], time_interval_range[1] >= 0
  )

  # Normalize actors_per_group to range
  if (length(actors_per_group) == 1) {
    actors_per_group <- c(actors_per_group, actors_per_group)
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

  # --- Determine Actions ---
  if (is.null(actions)) {
    if ("group_regulation" %in% action_categories) {
      # Built-in group regulation actions (SSRL)
      group_reg_actions <- c(
        "adapt", "cohesion", "consensus", "coregulate",
        "discuss", "emotion", "monitor", "plan", "synthesis"
      )
      if (length(action_categories) == 1 && action_categories == "group_regulation") {
        actions <- group_reg_actions
      } else {
        # Combine with other categories
        other_cats <- setdiff(action_categories, "group_regulation")
        other_actions <- if (length(other_cats) > 0) {
          get_learning_states(other_cats, n = max(1, n_actions - length(group_reg_actions)))
        } else {
          character(0)
        }
        actions <- unique(c(group_reg_actions, other_actions))
      }
    } else {
      # Use learning states from specified categories
      actions <- get_learning_states(action_categories, n = n_actions)
    }
  }

  n_action_types <- length(actions)
  stopifnot(n_action_types >= 2)

  # --- Generate Transition Probabilities ---
  if (is.null(transition_probs)) {
    transition_probs <- seqHMM::simulate_transition_probs(n_action_types)
  }
  dimnames(transition_probs) <- list(actions, actions)

  if (is.null(initial_probs)) {
    initial_probs <- seqHMM::simulate_initial_probs(n_action_types)
  }
  names(initial_probs) <- actions

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
    if (actors_per_group[1] == actors_per_group[2]) {
      n_actors_in_group <- actors_per_group[1]
    } else {
      n_actors_in_group <- sample(actors_per_group[1]:actors_per_group[2], 1)
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
      sequence[1] <- sample(actions, 1, prob = initial_probs)

      if (seq_len > 1) {
        for (j in 2:seq_len) {
          current_state <- sequence[j - 1]
          current_idx <- match(current_state, actions)
          sequence[j] <- sample(actions, 1, prob = transition_probs[current_idx, ])
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
#'   actions = GROUP_REGULATION_ACTIONS,
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
