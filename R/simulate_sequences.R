#' Simulate Markov Chain Sequences (Basic)
#'
#' @description
#' Generate sequences of states from a Markov chain. Can either use provided
#' transition matrix and initial probabilities, or auto-generate random ones
#' with optional learning state names.
#'
#' @param transition_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names. If NULL, random
#'   probabilities are generated based on `n_states`. Default: NULL.
#' @param initial_probabilities Named numeric vector of initial state
#'   probabilities. Must sum to 1. If NULL, random probabilities are generated.
#'   Default: NULL.
#' @param max_seq_length Integer. Maximum length of each sequence. Default: 20.
#' @param num_rows Integer. Number of sequences (rows) to generate. Default: 100.
#' @param n_states Integer. Number of states when auto-generating probabilities.
#'   Ignored if `transition_matrix` is provided. Default: 5.
#' @param state_names Character vector. Names for states when auto-generating.
#'   If NULL, uses letters (A, B, C, ...) or learning states if enabled.
#'   Ignored if `transition_matrix` is provided. Default: NULL.
#' @param use_learning_states Logical. If TRUE and auto-generating, uses
#'   realistic learning action verbs as state names. Default: FALSE.
#' @param learning_categories Character vector. Categories of learning states
#'   to use. Options: "metacognitive", "cognitive", "behavioral", "social",
#'   "motivational", "affective", "group_regulation", or "all".
#'   Only used if `use_learning_states = TRUE`. Default: "all".
#' @param smart_select Logical. If TRUE, intelligently selects learning states
#'   based on n_states. Only used if `use_learning_states = TRUE`. Default: TRUE.
#' @param min_na Integer. Minimum number of NA values per sequence. Default: 0.
#' @param max_na Integer. Maximum number of NA values per sequence. Default: 0.
#' @param include_na Logical. Whether to include NAs in sequences. Default: TRUE.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#' @param return_params Logical. If TRUE, returns a list with sequences and
#'   the generating parameters. Default: FALSE.
#'
#' @return If `return_params = FALSE` (default), a data frame with `num_rows`
#'   rows and `max_seq_length` columns. Each row is a sequence of state names,
#'   potentially with trailing NAs.
#'
#'   If `return_params = TRUE`, a list containing:
#'   \itemize{
#'     \item sequences: The data frame of sequences.
#'     \item transition_matrix: The transition probability matrix used.
#'     \item initial_probabilities: The initial probabilities used.
#'     \item state_names: The state names used.
#'   }
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
#' **Auto-generation mode**: When `transition_matrix` is NULL, the function
#' automatically generates random transition probabilities and initial
#' probabilities using `seqHMM::simulate_transition_probs()` and
#' `seqHMM::simulate_initial_probs()`.
#'
#' **Learning States**: When `use_learning_states = TRUE`, state names are
#' drawn from a curated collection of 180+ student learning action verbs:
#' \itemize{
#'   \item **metacognitive**: Plan, Monitor, Evaluate, Reflect, ...
#'   \item **cognitive**: Read, Study, Analyze, Summarize, ...
#'   \item **behavioral**: Practice, Annotate, Research, Review, ...
#'   \item **social**: Collaborate, Discuss, Explain, Share, ...
#'   \item **motivational**: Focus, Persist, Explore, Strive, ...
#'   \item **affective**: Enjoy, Appreciate, Cope, Curious, ...
#'   \item **group_regulation**: Adapt, Cohesion, Consensus, ...
#' }
#'
#' @examples
#' \dontrun{
#' # Simplest usage: all defaults (100 sequences, 20 length, 5 states)
#' sequences <- simulate_sequences(seed = 42)
#'
#' # With learning states
#' sequences <- simulate_sequences(use_learning_states = TRUE, seed = 42)
#'
#' # Method 1: Provide your own transition matrix
#' trans_mat <- matrix(c(
#'   0.7, 0.2, 0.1,
#'   0.3, 0.5, 0.2,
#'   0.2, 0.3, 0.5
#' ), nrow = 3, byrow = TRUE)
#' rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
#' init_probs <- c(A = 0.5, B = 0.3, C = 0.2)
#'
#' sequences <- simulate_sequences(
#'   transition_matrix = trans_mat,
#'   initial_probabilities = init_probs,
#'   max_seq_length = 20,
#'   num_rows = 100
#' )
#'
#' # Method 2: Auto-generate with letter names
#' sequences <- simulate_sequences(
#'   max_seq_length = 20,
#'   num_rows = 100,
#'   n_states = 5,
#'   seed = 42
#' )
#'
#' # Method 3: Auto-generate with learning state names
#' sequences <- simulate_sequences(
#'   max_seq_length = 25,
#'   num_rows = 150,
#'   n_states = 6,
#'   use_learning_states = TRUE,
#'   learning_categories = c("metacognitive", "cognitive"),
#'   seed = 123
#' )
#'
#' # Method 4: Get sequences AND the generating parameters
#' result <- simulate_sequences(
#'   max_seq_length = 20,
#'   num_rows = 100,
#'   n_states = 4,
#'   use_learning_states = TRUE,
#'   return_params = TRUE,
#'   seed = 42
#' )
#' result$sequences          # The sequences
#' result$transition_matrix  # The random transition matrix
#' result$state_names        # e.g., c("Plan", "Monitor", "Read", "Practice")
#'
#' # Method 5: Custom state names
#' sequences <- simulate_sequences(
#'   max_seq_length = 20,
#'   num_rows = 100,
#'   n_states = 4,
#'   state_names = c("Explore", "Learn", "Practice", "Master"),
#'   seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_sequences_advanced}} for sequences with stability modes,
#' \code{\link{simulate_long_data}} for long-format educational data,
#' \code{\link{get_learning_states}} for available learning state verbs,
#' \code{\link{generate_tna_networks}} for generating complete TNA networks.
#'
#' @importFrom stats runif
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @export
simulate_sequences <- function(transition_matrix = NULL,
                                initial_probabilities = NULL,
                                max_seq_length = 20,
                                num_rows = 100,
                                n_states = 5,
                                state_names = NULL,
                                use_learning_states = FALSE,
                                learning_categories = "all",
                                smart_select = TRUE,
                                min_na = 0,
                                max_na = 0,
                                include_na = TRUE,
                                seed = NULL,
                                return_params = FALSE) {
  # Set seed if provided

if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Auto-generate probabilities if not provided ---
  auto_generated <- FALSE
  if (is.null(transition_matrix) || is.null(initial_probabilities)) {
    auto_generated <- TRUE

    # Determine state names
    if (use_learning_states) {
      if (smart_select && ("all" %in% learning_categories || length(learning_categories) > 2)) {
        state_names <- smart_select_states(
          n_states = n_states,
          primary_categories = if ("all" %in% learning_categories) NULL else learning_categories
        )
      } else {
        state_names <- get_learning_states(
          categories = learning_categories,
          n = n_states
        )
      }
    } else if (is.null(state_names)) {
      # Default to letters or S1, S2, ...
      if (n_states <= 26) {
        state_names <- LETTERS[1:n_states]
      } else {
        state_names <- paste0("S", 1:n_states)
      }
    } else {
      # Use provided state_names
      if (length(state_names) < n_states) {
        stop("state_names must have at least n_states elements.")
      }
      state_names <- state_names[1:n_states]
    }

    # Generate random probabilities
    if (is.null(transition_matrix)) {
      transition_matrix <- seqHMM::simulate_transition_probs(n_states)
      dimnames(transition_matrix) <- list(state_names, state_names)
    }

    if (is.null(initial_probabilities)) {
      initial_probabilities <- seqHMM::simulate_initial_probs(n_states)
      names(initial_probabilities) <- state_names
    }
  }

  # --- Input validation ---
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
  num_states_actual <- length(states_names)

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
    current_state_index <- sample(num_states_actual, 1, prob = initial_probabilities)
    sequence[1] <- states_names[current_state_index]

    # Generate remaining states in sequence
    if (seq_length > 1) {
      for (i in 2:seq_length) {
        current_state_index <- match(sequence[i - 1], states_names)
        next_state_index <- sample(num_states_actual, 1, prob = transition_matrix[current_state_index, ])
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

  # Return based on return_params
  if (return_params) {
    return(list(
      sequences = result_df,
      transition_matrix = transition_matrix,
      initial_probabilities = initial_probabilities,
      state_names = states_names
    ))
  }

  return(result_df)
}
