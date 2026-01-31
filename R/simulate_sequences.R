#' Simulate Markov Chain Sequences (Basic)
#'
#' @description
#' Generate sequences of states from a Markov chain given a transition matrix
#' and initial state probabilities. This is the basic version without
#' stability modes.
#'
#' @param transition_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names.
#' @param initial_probabilities Named numeric vector of initial state
#'   probabilities. Must sum to 1.
#' @param max_seq_length Integer. Maximum length of each sequence.
#' @param num_rows Integer. Number of sequences (rows) to generate.
#' @param min_na Integer. Minimum number of NA values per sequence (default: 0).
#' @param max_na Integer. Maximum number of NA values per sequence (default: 0).
#' @param include_na Logical. Whether to include NAs in sequences (default: TRUE).
#'
#' @return A data frame with `num_rows` rows and `max_seq_length` columns.
#'   Each row is a sequence of state names, potentially with trailing NAs.
#'
#' @details
#' The function generates sequences by:
#' \enumerate{
#'   \item Sampling an initial state from `initial_probabilities`.
#'   \item For each subsequent position, sampling the next state based on
#'     the current state's row in `transition_matrix`.
#'   \item If `include_na` is TRUE, trailing NAs are added to achieve a
#'     random number of NAs between `min_na` and `max_na`.
#' }
#'
#' Input validation ensures:
#' \itemize{
#'   \item Transition matrix rows sum to 1 (tolerance: 1e-6).
#'   \item Initial probabilities sum to 1 (tolerance: 1e-6).
#'   \item Valid NA range parameters.
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple 3-state transition matrix
#' trans_mat <- matrix(c(
#'   0.7, 0.2, 0.1,
#'   0.3, 0.5, 0.2,
#'   0.2, 0.3, 0.5
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#'
#' init_probs <- c(A = 0.5, B = 0.3, C = 0.2)
#'
#' # Generate 100 sequences of length 20
#' sequences <- simulate_sequences(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 20,
#'   num_rows = 100,
#'   min_na = 0,
#'   max_na = 5
#' )
#' }
#'
#' @importFrom stats runif
#' @export
simulate_sequences <- function(transition_matrix,
                                initial_probabilities,
                                max_seq_length,
                                num_rows,
                                min_na = 0,
                                max_na = 0,
                                include_na = TRUE) {
  # Input validation
  if (is.null(transition_matrix) || is.null(initial_probabilities)) {
    stop("Transition matrix or initial probabilities are not defined.")
  }

  # Ensure the transition matrix rows sum to 1
  if (!all(abs(rowSums(transition_matrix) - 1) < 1e-6)) {
    stop("Rows of the transition matrix must sum to 1.")
  }

  # Ensure the initial probabilities sum to 1
  if (abs(sum(initial_probabilities) - 1) >= 1e-6) {
    stop("Initial probabilities must sum to 1.")
  }

  # Check if min_na and max_na are valid
  if (min_na < 0 || max_na > max_seq_length || min_na > max_na) {
    stop("Invalid 'min_na' or 'max_na' values. 'min_na' must be >=0, 'max_na' must be <= max_seq_length and 'min_na' must be <= 'max_na'.")
  }

  states_names <- rownames(transition_matrix)
  num_states <- length(states_names)

  # Pre-allocate result matrix
  result_matrix <- matrix(NA, nrow = num_rows, ncol = max_seq_length)

  for (row in 1:num_rows) {
    if (include_na) {
      # Determine number of NAs for this sequence
      num_na <- if (min_na == max_na) min_na else sample(min_na:max_na, 1)

      # Calculate actual sequence length (non-NA values)
      seq_length <- max_seq_length - num_na
    } else {
      seq_length <- max_seq_length
    }

    # Ensure sequence length is at least 1 if the entire row shouldn't be NA
    seq_length <- max(1, seq_length)

    # Generate complete sequence first
    sequence <- character(seq_length)

    # Sample initial state
    current_state_index <- sample(num_states, 1, prob = initial_probabilities)
    sequence[1] <- states_names[current_state_index]

    # Generate remaining states in sequence
    if (seq_length > 1) {
      for (i in 2:seq_length) {
        current_state_index <- match(sequence[i - 1], states_names)
        next_state_index <- sample(num_states, 1, prob = transition_matrix[current_state_index, ])
        sequence[i] <- states_names[next_state_index]
      }
    }

    # Fill the row: first with valid states, then with NAs
    result_matrix[row, 1:seq_length] <- sequence
    # The rest of the row will remain NA (initialized that way)
  }

  # Convert to data frame with standard column names
  result_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
  names(result_df) <- paste0("V", 1:max_seq_length)

  return(result_df)
}
