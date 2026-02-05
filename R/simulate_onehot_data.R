#' Simulate One-Hot Encoded Sequence Data
#'
#' @description
#' Generate simulated sequence data in one-hot encoded format. Creates the same
#' hierarchical structure as \code{\link{simulate_long_data}} (Actor, Achiever,
#' Group, Course, Time) but with binary 0/1 columns for each state instead of
#' a categorical Action column.
#'
#' @inheritParams simulate_long_data
#' @param sort_states Logical. Sort state columns alphabetically. Default: FALSE.
#' @param state_prefix Character. Prefix for state column names. Default: "".
#'
#' @return A tibble (data frame) with columns:
#' \describe{
#'   \item{Actor}{Integer. Actor identifier (1 to n_actors).}
#'   \item{Achiever}{Character. Achievement level (e.g., "High", "Low").}
#'   \item{Group}{Numeric. Group identifier (1 to n_groups).}
#'   \item{Course}{Character. Course identifier.}
#'   \item{Time}{POSIXct. Timestamp of the action.}
#'   \item{<state columns>}{Integer (0/1). One column per state, with 1 indicating
#'     the active state for that row.}
#' }
#'
#' @details
#' This function internally calls \code{\link{simulate_long_data}} to generate
#' the base data, then applies \code{\link{action_to_onehot}} to convert the
#' Action column to one-hot encoded columns.
#'
#' Each row will have exactly one state column with value 1, and all others
#' with value 0, representing the action performed at that time point.
#'
#' @examples
#' \dontrun{
#' # Basic usage with new standardized params
#' data <- simulate_onehot_data(n_groups = 5, n_actors = 3, seed = 42)
#' head(data)
#'
#' # Custom states with prefix
#' data <- simulate_onehot_data(
#'   n_groups = 3,
#'   states = c("Read", "Write", "Think"),
#'   state_prefix = "state_",
#'   seed = 123
#' )
#' names(data)
#'
#' # Verify one-hot encoding (each row sums to 1)
#' state_cols <- setdiff(names(data), c("Actor", "Achiever", "Group", "Course", "Time"))
#' all(rowSums(data[, state_cols]) == 1)  # TRUE
#'
#' # Old parameter names still work (backward compatible)
#' data <- simulate_onehot_data(
#'   actors_per_group = 5,
#'   actions = c("A", "B", "C"),
#'   seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_long_data}} for categorical Action format,
#' \code{\link{action_to_onehot}} for converting existing long data,
#' \code{\link{simulate_sequences}} for generating wide-format sequences.
#'
#' @export
simulate_onehot_data <- function(
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
    # Backward compatibility - old parameter names
    actors_per_group = NULL,
    actions = NULL,
    action_categories = NULL,
    n_actions = NULL,
    transition_probs = NULL,
    initial_probs = NULL
) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(actors_per_group)) n_actors <- actors_per_group
  if (!is.null(actions)) states <- actions
  if (!is.null(action_categories)) categories <- action_categories
  if (!is.null(n_actions)) n_states <- n_actions
  if (!is.null(transition_probs)) trans_matrix <- transition_probs
  if (!is.null(initial_probs)) init_probs <- initial_probs

  # Generate long format data
  long_data <- simulate_long_data(
    n_groups = n_groups,
    n_actors = n_actors,
    n_courses = n_courses,
    n_states = n_states,
    states = states,
    use_learning_states = use_learning_states,
    categories = categories,
    seq_length_range = seq_length_range,
    achiever_levels = achiever_levels,
    achiever_probs = achiever_probs,
    start_time = start_time,
    time_interval_range = time_interval_range,
    trans_matrix = trans_matrix,
    init_probs = init_probs,
    seed = seed
  )

  # Convert to one-hot encoding
  result <- action_to_onehot(
    data = long_data,
    action_col = "Action",
    states = NULL,
    drop_action = TRUE,
    sort_states = sort_states,
    prefix = state_prefix
  )

  # Reorder columns: metadata first, then state columns
  metadata_cols <- c("Actor", "Achiever", "Group", "Course", "Time")
  state_cols <- setdiff(names(result), metadata_cols)
  result <- result[, c(metadata_cols, state_cols), drop = FALSE]

  # Preserve tibble class
  result <- structure(
    result,
    class = c("tbl_df", "tbl", "data.frame")
  )

  return(result)
}
