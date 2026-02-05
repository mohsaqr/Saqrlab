#' Simulate Markov Chain Sequences (Advanced)
#'
#' @description
#' Generate sequences of states from a Markov chain with advanced stability
#' and instability modes. Supports stable transitions, probability perturbation,
#' and unlikely jump modes.
#'
#' @param n_sequences Integer. Number of sequences (rows) to generate. Default: 1000.
#' @param seq_length Integer. Maximum length of each sequence. Default: 20.
#' @param n_states Integer. Number of states when auto-generating probabilities.
#'   Ignored if `trans_matrix` is provided. Default: 5.
#' @param states Character vector. Names for states when auto-generating.
#'   If NULL, uses letters (A, B, C, ...) or learning states if enabled.
#'   Ignored if `trans_matrix` is provided. Default: NULL.
#' @param use_learning_states Logical. If TRUE and auto-generating, uses
#'   realistic learning action verbs as state names. Default: TRUE.
#' @param categories Character vector. Categories of learning states to use.
#'   Only used if `use_learning_states = TRUE`. Default: "all".
#' @param trans_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names. If NULL, random
#'   probabilities are generated. Default: NULL.
#' @param init_probs Named numeric vector of initial state probabilities.
#'   Must sum to 1. If NULL, random probabilities are generated. Default: NULL.
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
#' @param na_range Integer vector of length 2 (min, max) or single integer (min=max).
#'   Range of NA values per sequence. Default: c(0, 0).
#' @param include_na Logical. Whether to include NAs in sequences. Default: TRUE.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @param transition_matrix Deprecated. Use `trans_matrix` instead.
#' @param initial_probabilities Deprecated. Use `init_probs` instead.
#' @param num_rows Deprecated. Use `n_sequences` instead.
#' @param max_seq_length Deprecated. Use `seq_length` instead.
#' @param min_na Deprecated. Use `na_range` instead.
#' @param max_na Deprecated. Use `na_range` instead.
#'
#' @return A data frame with `n_sequences` rows and `seq_length` columns.
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
#' # Simplest usage: all defaults
#' sequences <- simulate_sequences_advanced(seed = 42)
#'
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
#'   trans_matrix = trans_mat,
#'   init_probs = init_probs,
#'   seq_length = 30,
#'   n_sequences = 100,
#'   stable_transitions = stable,
#'   stability_prob = 0.95,
#'   unstable_mode = "unlikely_jump",
#'   unstable_random_transition_prob = 0.3,
#'   na_range = c(0, 5)
#' )
#'
#' # Old parameter names still work (backward compatible)
#' sequences <- simulate_sequences_advanced(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 30,
#'   num_rows = 100
#' )
#' }
#'
#' @importFrom stats runif
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @export
simulate_sequences_advanced <- function(n_sequences = 1000,
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
                                         # Backward compatibility - old parameter names
                                         transition_matrix = NULL,
                                         initial_probabilities = NULL,
                                         num_rows = NULL,
                                         max_seq_length = NULL,
                                         min_na = NULL,
                                         max_na = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(transition_matrix)) trans_matrix <- transition_matrix
  if (!is.null(initial_probabilities)) init_probs <- initial_probabilities
  if (!is.null(num_rows)) n_sequences <- num_rows
  if (!is.null(max_seq_length)) seq_length <- max_seq_length

  # Handle na_range from min_na/max_na
  if (!is.null(min_na) || !is.null(max_na)) {
    min_na_val <- if (!is.null(min_na)) min_na else 0
    max_na_val <- if (!is.null(max_na)) max_na else min_na_val
    na_range <- c(min_na_val, max_na_val)
  }

  # Normalize na_range to length 2
  if (length(na_range) == 1) {
    na_range <- c(na_range, na_range)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Auto-generate probabilities if not provided ---
  if (is.null(trans_matrix) || is.null(init_probs)) {
    # Determine state names
    if (use_learning_states) {
      states <- get_learning_states(categories = categories, n = n_states)
    } else if (is.null(states)) {
      if (n_states <= 26) {
        states <- LETTERS[1:n_states]
      } else {
        states <- paste0("S", 1:n_states)
      }
    } else {
      if (length(states) < n_states) {
        stop("states must have at least n_states elements.")
      }
      states <- states[1:n_states]
    }

    # Generate random probabilities
    if (is.null(trans_matrix)) {
      trans_matrix <- seqHMM::simulate_transition_probs(n_states)
      dimnames(trans_matrix) <- list(states, states)
    }

    if (is.null(init_probs)) {
      init_probs <- seqHMM::simulate_initial_probs(n_states)
      names(init_probs) <- states
    }
  }

  # --- Input Validation ---
  stopifnot(
    is.matrix(trans_matrix), is.numeric(trans_matrix), all(!is.na(trans_matrix)),
    nrow(trans_matrix) == ncol(trans_matrix),
    is.numeric(init_probs), length(init_probs) == nrow(trans_matrix),
    abs(sum(init_probs) - 1) < 1e-6, all(init_probs >= 0),
    all(abs(rowSums(trans_matrix) - 1) < 1e-6), all(trans_matrix >= 0),
    is.numeric(seq_length), length(seq_length) == 1, seq_length > 0,
    is.numeric(n_sequences), length(n_sequences) == 1, n_sequences > 0,
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
    is.numeric(na_range), length(na_range) == 2,
    na_range[1] >= 0, na_range[2] >= 0, na_range[1] <= na_range[2],
    is.logical(include_na), length(include_na) == 1
  )

  states_names <- rownames(trans_matrix)
  if (is.null(states_names) || any(is.na(states_names)) || anyDuplicated(states_names)) {
    stop("`trans_matrix` must have unique, non-NA row names.")
  }
  num_states <- length(states_names)

  # Cap max_na
  max_allowable_na <- max(0, seq_length - 2)
  effective_max_na <- min(na_range[2], max_allowable_na)
  effective_min_na <- min(na_range[1], effective_max_na)

  if (effective_min_na > effective_max_na && na_range[2] > 0) {
    warning(sprintf(
      "NA range invalid (min=%d > max=%d) for seq_length=%d; setting NAs to 0.",
      effective_min_na, effective_max_na, seq_length
    ))
    effective_min_na <- effective_max_na <- 0
  }

  # Pre-allocate & Pre-process
  result_matrix <- matrix(NA_character_, nrow = n_sequences, ncol = seq_length)
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
  result_matrix[, 1] <- sample(states_names, n_sequences, replace = TRUE, prob = init_probs)
  for (j in 2:seq_length) {
    for (i in 1:n_sequences) {
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
        original_probs <- trans_matrix[current_state_index, ]
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
  if (include_na && effective_max_na > 0 && seq_length >= 2) {
    for (row in 1:n_sequences) {
      num_na <- if (effective_min_na == effective_max_na) effective_min_na else sample(effective_min_na:effective_max_na, 1)
      if (num_na > 0) {
        na_positions <- tail(1:seq_length, num_na)
        result_matrix[row, na_positions] <- NA_character_
      }
    }
  } else if (include_na && na_range[2] > 0 && seq_length < 2) {
    warning("Cannot add NAs while preserving 2 elements when seq_length < 2.")
  }

  result_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
  names(result_df) <- paste0("V", 1:seq_length)
  return(result_df)
}
