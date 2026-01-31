#' Simulate Markov Chain Sequences (Advanced)
#'
#' @description
#' Generate sequences of states from a Markov chain with advanced stability
#' and instability modes. Supports stable transitions, probability perturbation,
#' and unlikely jump modes.
#'
#' @param transition_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names.
#' @param initial_probabilities Named numeric vector of initial state
#'   probabilities. Must sum to 1.
#' @param max_seq_length Integer. Maximum length of each sequence.
#' @param num_rows Integer. Number of sequences (rows) to generate.
#' @param stable_transitions List of character vectors. Each vector contains
#'   two state names defining a stable transition pair (from, to).
#'   Default: NULL (no stable transitions).
#' @param stability_prob Numeric in (0 to 1). Probability of following a stable
#'   transition when in a stable state. Default: 0.95.
#' @param unstable_mode Character. Mode for unstable transitions. One of:
#'   \itemize{
#'     \item "random_jump": Uniform random jump to any state.
#'     \item "perturb_prob": Perturb transition probabilities with noise.
#'     \item "unlikely_jump": Jump to states with low probability.
#'   }
#'   Default: "unlikely_jump".
#' @param unstable_random_transition_prob Numeric in (0 to 1). Probability of
#'   taking an unstable action. Default: 0.4.
#' @param unstable_perturb_noise Numeric in (0 to 1). Noise factor for probability
#'   perturbation mode. Default: 0.5.
#' @param unlikely_prob_threshold Numeric in (0 to 1). Threshold below which
#'   transitions are considered "unlikely". Default: 0.1.
#' @param min_na Integer. Minimum number of NA values per sequence. Default: 0.
#' @param max_na Integer. Maximum number of NA values per sequence. Default: 0.
#' @param include_na Logical. Whether to include NAs in sequences. Default: TRUE.
#'
#' @return A data frame with `num_rows` rows and `max_seq_length` columns.
#'   Each row is a sequence of state names, potentially with trailing NAs.
#'
#' @details
#' The function extends basic Markov chain simulation with:
#'
#' **Stable Transitions**: If the current state has a defined stable transition
#' and a random draw is below `stability_prob`, the sequence follows the stable
#' transition directly.
#'
#' **Unstable Modes** (when not following stable transitions):
#' \itemize{
#'   \item "random_jump": With probability `unstable_random_transition_prob`,
#'     jump uniformly to any state.
#'   \item "perturb_prob": Multiply transition probabilities by random noise
#'     in range `[1-noise, 1+noise]`, then renormalize.
#'   \item "unlikely_jump": With probability `unstable_random_transition_prob`,
#'     jump to a state with transition probability below threshold.
#' }
#'
#' NAs are added at the end of sequences, preserving at least 2 non-NA values.
#'
#' @examples
#' \dontrun{
#' # Create a 4-state transition matrix
#' trans_mat <- matrix(c(
#'   0.6, 0.2, 0.1, 0.1,
#'   0.1, 0.7, 0.1, 0.1,
#'   0.1, 0.1, 0.6, 0.2,
#'   0.2, 0.1, 0.1, 0.6
#' ), nrow = 4, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C", "D")
#'
#' init_probs <- c(A = 0.25, B = 0.25, C = 0.25, D = 0.25)
#'
#' # Define stable transitions: A->B and C->D are "stable"
#' stable <- list(c("A", "B"), c("C", "D"))
#'
#' # Generate sequences with unlikely_jump mode
#' sequences <- simulate_sequences_advanced(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 30,
#'   num_rows = 100,
#'   stable_transitions = stable,
#'   stability_prob = 0.95,
#'   unstable_mode = "unlikely_jump",
#'   unstable_random_transition_prob = 0.3,
#'   min_na = 0,
#'   max_na = 5
#' )
#' }
#'
#' @importFrom stats runif
#' @export
simulate_sequences_advanced <- function(transition_matrix,
                                         initial_probabilities,
                                         max_seq_length,
                                         num_rows,
                                         stable_transitions = NULL,
                                         stability_prob = 0.95,
                                         unstable_mode = "unlikely_jump",
                                         unstable_random_transition_prob = 0.4,
                                         unstable_perturb_noise = 0.5,
                                         unlikely_prob_threshold = 0.1,
                                         min_na = 0,
                                         max_na = 0,
                                         include_na = TRUE) {
  # --- Input Validation ---
  stopifnot(
    is.matrix(transition_matrix), is.numeric(transition_matrix), all(!is.na(transition_matrix)),
    nrow(transition_matrix) == ncol(transition_matrix),
    is.numeric(initial_probabilities), length(initial_probabilities) == nrow(transition_matrix),
    abs(sum(initial_probabilities) - 1) < 1e-6, all(initial_probabilities >= 0),
    all(abs(rowSums(transition_matrix) - 1) < 1e-6), all(transition_matrix >= 0),
    is.numeric(max_seq_length), length(max_seq_length) == 1, max_seq_length > 0,
    is.numeric(num_rows), length(num_rows) == 1, num_rows > 0,
    is.null(stable_transitions) || is.list(stable_transitions),
    is.numeric(stability_prob), length(stability_prob) == 1, stability_prob >= 0, stability_prob <= 1,
    is.character(unstable_mode), length(unstable_mode) == 1,
    unstable_mode %in% c("random_jump", "perturb_prob", "unlikely_jump"),
    is.numeric(unstable_random_transition_prob), length(unstable_random_transition_prob) == 1,
    unstable_random_transition_prob >= 0, unstable_random_transition_prob <= 1,
    is.numeric(unstable_perturb_noise), length(unstable_perturb_noise) == 1,
    unstable_perturb_noise >= 0, unstable_perturb_noise <= 1,
    is.numeric(unlikely_prob_threshold), length(unlikely_prob_threshold) == 1,
    unlikely_prob_threshold >= 0, unlikely_prob_threshold <= 1,
    is.numeric(min_na), length(min_na) == 1, min_na >= 0,
    is.numeric(max_na), length(max_na) == 1, max_na >= 0, min_na <= max_na,
    is.logical(include_na), length(include_na) == 1
  )

  states_names <- rownames(transition_matrix)
  if (is.null(states_names) || any(is.na(states_names)) || anyDuplicated(states_names)) {
    stop("`transition_matrix` must have unique, non-NA row names.")
  }
  num_states <- length(states_names)

  # Cap max_na
  max_allowable_na <- max(0, max_seq_length - 2)
  effective_max_na <- min(max_na, max_allowable_na)
  effective_min_na <- min(min_na, effective_max_na)

  if (effective_min_na > effective_max_na && max_na > 0) {
    warning(sprintf(
      "NA range invalid (min=%d > max=%d) for max_seq_length=%d; setting NAs to 0.",
      effective_min_na, effective_max_na, max_seq_length
    ))
    effective_min_na <- effective_max_na <- 0
  }

  # Pre-allocate & Pre-process
  result_matrix <- matrix(NA_character_, nrow = num_rows, ncol = max_seq_length)
  stable_lookup <- list()
  if (!is.null(stable_transitions)) {
    for (pair in stable_transitions) {
      if (length(pair) == 2 && all(pair %in% states_names)) {
        stable_lookup[[pair[1]]] <- pair[2]
      } else {
        warning(paste("Invalid stable transition pair:", paste(pair, collapse = "->"), "- skipping"))
      }
    }
  }

  # --- Generate Sequences ---
  result_matrix[, 1] <- sample(states_names, num_rows, replace = TRUE, prob = initial_probabilities)
  for (j in 2:max_seq_length) {
    for (i in 1:num_rows) {
      current_state <- result_matrix[i, j - 1]
      if (is.na(current_state)) {
        result_matrix[i, j] <- NA_character_
        next
      }
      current_state_index <- match(current_state, states_names)
      if (is.na(current_state_index)) {
        warning(sprintf(
          "State '%s' not found at row %d, col %d. Assigning NA.",
          current_state, i, j
        ))
        result_matrix[i, j] <- NA_character_
        next
      }

      if (current_state %in% names(stable_lookup) && runif(1) < stability_prob) {
        result_matrix[i, j] <- stable_lookup[[current_state]]
      } else {
        original_probs <- transition_matrix[current_state_index, ]
        if (unstable_mode == "perturb_prob") {
          noise <- runif(num_states, 1 - unstable_perturb_noise, 1 + unstable_perturb_noise)
          perturbed_probs <- pmax(0, original_probs * noise)
          sum_probs <- sum(perturbed_probs)
          final_probs <- if (sum_probs > 1e-9) perturbed_probs / sum_probs else rep(1 / num_states, num_states)
          result_matrix[i, j] <- sample(states_names, 1, prob = final_probs)
        } else if (unstable_mode == "unlikely_jump") {
          if (runif(1) < unstable_random_transition_prob) {
            unlikely_indices <- which(original_probs < unlikely_prob_threshold)
            result_matrix[i, j] <- if (length(unlikely_indices) > 0) states_names[sample(unlikely_indices, 1)] else sample(states_names, 1)
          } else {
            result_matrix[i, j] <- sample(states_names, 1, prob = original_probs)
          }
        } else {
          # Default: "random_jump"
          if (runif(1) < unstable_random_transition_prob) {
            result_matrix[i, j] <- sample(states_names, 1) # Uniform jump
          } else {
            result_matrix[i, j] <- sample(states_names, 1, prob = original_probs) # Original probs
          }
        }
      }
    }
  }

  # --- Post-processing for NAs ---
  if (include_na && effective_max_na > 0 && max_seq_length >= 2) {
    for (row in 1:num_rows) {
      num_na <- if (effective_min_na == effective_max_na) effective_min_na else sample(effective_min_na:effective_max_na, 1)
      if (num_na > 0) {
        na_positions <- tail(1:max_seq_length, num_na)
        result_matrix[row, na_positions] <- NA_character_
      }
    }
  } else if (include_na && max_na > 0 && max_seq_length < 2) {
    warning("Cannot add NAs while preserving 2 elements when max_seq_length < 2.")
  }

  result_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
  names(result_df) <- paste0("V", 1:max_seq_length)
  return(result_df)
}
