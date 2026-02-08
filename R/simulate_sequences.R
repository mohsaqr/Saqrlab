#' Simulate Markov Chain Sequences (Basic)
#'
#' @description
#' Generate sequences of states from a Markov chain. Can either use provided
#' transition matrix and initial probabilities, or auto-generate random ones
#' with optional learning state names.
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
#' @param categories Character vector. Categories of learning states
#'   to use. Options: "metacognitive", "cognitive", "behavioral", "social",
#'   "motivational", "affective", "group_regulation", or "all".
#'   Only used if `use_learning_states = TRUE`. Default: "all".
#' @param smart_select Logical. If TRUE, intelligently selects learning states
#'   based on n_states. Only used if `use_learning_states = TRUE`. Default: TRUE.
#' @param trans_matrix Square numeric matrix of transition probabilities.
#'   Rows must sum to 1. Row names define state names. If NULL, random
#'   probabilities are generated based on `n_states`. Default: NULL.
#' @param init_probs Named numeric vector of initial state
#'   probabilities. Must sum to 1. If NULL, random probabilities are generated.
#'   Default: NULL.
#' @param na_range Integer vector of length 2 (min, max) or single integer (min=max).
#'   Range of NA values per sequence. Default: c(0, 0).
#' @param include_na Logical. Whether to include NAs in sequences. Default: TRUE.
#' @param include_params Logical. If TRUE, returns a list with sequences and
#'   the generating parameters. Default: FALSE.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @param num_rows Deprecated. Use `n_sequences` instead.
#' @param max_seq_length Deprecated. Use `seq_length` instead.
#' @param state_names Deprecated. Use `states` instead.
#' @param learning_categories Deprecated. Use `categories` instead.
#' @param transition_matrix Deprecated. Use `trans_matrix` instead.
#' @param initial_probabilities Deprecated. Use `init_probs` instead.
#' @param min_na Deprecated. Use `na_range` instead.
#' @param max_na Deprecated. Use `na_range` instead.
#' @param return_params Deprecated. Use `include_params` instead.
#'
#' @return If `include_params = FALSE` (default), a data frame with `n_sequences`
#'   rows and `seq_length` columns. Each row is a sequence of state names,
#'   potentially with trailing NAs.
#'
#'   If `include_params = TRUE`, a list containing:
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
#'   \item Sampling an initial state from `init_probs`.
#'   \item For each subsequent position, sampling the next state based on
#'     the current state's row in `trans_matrix`.
#'   \item If `include_na` is TRUE, trailing NAs are added to achieve a
#'     random number of NAs within `na_range`.
#' }
#'
#' **Auto-generation mode**: When `trans_matrix` is NULL, the function
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
#' # Simplest usage: all defaults (1000 sequences, 20 length, 5 states, learning states)
#' sequences <- simulate_sequences(seed = 42)
#'
#' # Explicit new parameter names
#' sequences <- simulate_sequences(
#'   n_sequences = 100,
#'   seq_length = 15,
#'   n_states = 6,
#'   use_learning_states = TRUE,
#'   categories = c("metacognitive", "cognitive"),
#'   seed = 42
#' )
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
#'   trans_matrix = trans_mat,
#'   init_probs = init_probs,
#'   seq_length = 20,
#'   n_sequences = 100
#' )
#'
#' # Method 2: Auto-generate with letter names
#' sequences <- simulate_sequences(
#'   seq_length = 20,
#'   n_sequences = 100,
#'   n_states = 5,
#'   use_learning_states = FALSE,
#'   seed = 42
#' )
#'
#' # Method 3: Auto-generate with learning state names
#' sequences <- simulate_sequences(
#'   seq_length = 25,
#'   n_sequences = 150,
#'   n_states = 6,
#'   use_learning_states = TRUE,
#'   categories = c("metacognitive", "cognitive"),
#'   seed = 123
#' )
#'
#' # Method 4: Get sequences AND the generating parameters
#' result <- simulate_sequences(
#'   seq_length = 20,
#'   n_sequences = 100,
#'   n_states = 4,
#'   use_learning_states = TRUE,
#'   include_params = TRUE,
#'   seed = 42
#' )
#' result$sequences          # The sequences
#' result$transition_matrix  # The random transition matrix
#' result$state_names        # e.g., c("Plan", "Monitor", "Read", "Practice")
#'
#' # Method 5: Custom state names
#' sequences <- simulate_sequences(
#'   seq_length = 20,
#'   n_sequences = 100,
#'   n_states = 4,
#'   states = c("Explore", "Learn", "Practice", "Master"),
#'   seed = 42
#' )
#'
#' # Old parameter names still work (backward compatible)
#' sequences <- simulate_sequences(
#'   num_rows = 100,
#'   max_seq_length = 15,
#'   state_names = c("A", "B", "C"),
#'   seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_sequences_advanced}} for sequences with stability modes,
#' \code{\link{simulate_long_data}} for long-format educational data,
#' \code{\link{get_learning_states}} for available learning state verbs,
#' \code{\link{simulate_tna_networks}} for generating complete TNA networks.
#'
#' @importFrom stats runif
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @export
simulate_sequences <- function(n_sequences = 1000,
                                seq_length = 20,
                                n_states = 8,
                                states = NULL,
                                use_learning_states = TRUE,
                                categories = "all",
                                smart_select = TRUE,
                                trans_matrix = NULL,
                                init_probs = NULL,
                                na_range = c(0, 0),
                                include_na = TRUE,
                                include_params = FALSE,
                                seed = NULL,
                                # Backward compatibility - old parameter names
                                num_rows = NULL,
                                max_seq_length = NULL,
                                state_names = NULL,
                                learning_categories = NULL,
                                transition_matrix = NULL,
                                initial_probabilities = NULL,
                                min_na = NULL,
                                max_na = NULL,
                                return_params = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(num_rows)) n_sequences <- num_rows
  if (!is.null(max_seq_length)) seq_length <- max_seq_length
  if (!is.null(state_names)) states <- state_names
  if (!is.null(learning_categories)) categories <- learning_categories
  if (!is.null(transition_matrix)) trans_matrix <- transition_matrix
  if (!is.null(initial_probabilities)) init_probs <- initial_probabilities
  if (!is.null(return_params)) include_params <- return_params

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
  auto_generated <- FALSE
  if (is.null(trans_matrix) || is.null(init_probs)) {
    auto_generated <- TRUE

    # Determine state names
    if (use_learning_states) {
      if (smart_select && ("all" %in% categories || length(categories) > 2)) {
        states <- select_states(
          n_states = n_states,
          primary_categories = if ("all" %in% categories) NULL else categories
        )
      } else {
        states <- get_learning_states(
          categories = categories,
          n = n_states
        )
      }
    } else if (is.null(states)) {
      # Default to letters or S1, S2, ...
      if (n_states <= 26) {
        states <- LETTERS[1:n_states]
      } else {
        states <- paste0("S", 1:n_states)
      }
    } else {
      # Use provided states
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

  # --- Input validation ---
  # Ensure the transition matrix rows sum to 1
  if (!all(abs(rowSums(trans_matrix) - 1) < 1e-6)) {
    stop("Rows of the transition matrix must sum to 1.")
  }

  # Ensure the initial probabilities sum to 1
  if (abs(sum(init_probs) - 1) >= 1e-6) {
    stop("Initial probabilities must sum to 1.")
  }

  # Check if na_range is valid
  if (na_range[1] < 0 || na_range[2] > seq_length || na_range[1] > na_range[2]) {
    stop("Invalid 'na_range' values. min must be >=0, max must be <= seq_length, and min must be <= max.")
  }

  states_names <- rownames(trans_matrix)
  num_states_actual <- length(states_names)

  # Pre-allocate result matrix
  result_matrix <- matrix(NA, nrow = n_sequences, ncol = seq_length)

  for (row in 1:n_sequences) {
    if (include_na) {
      # Determine number of NAs for this sequence
      num_na <- if (na_range[1] == na_range[2]) na_range[1] else sample(na_range[1]:na_range[2], 1)

      # Calculate actual sequence length (non-NA values)
      actual_seq_length <- seq_length - num_na
    } else {
      actual_seq_length <- seq_length
    }

    # Ensure sequence length is at least 1 if the entire row shouldn't be NA
    actual_seq_length <- max(1, actual_seq_length)

    # Generate complete sequence first
    sequence <- character(actual_seq_length)

    # Sample initial state
    current_state_index <- sample(num_states_actual, 1, prob = init_probs)
    sequence[1] <- states_names[current_state_index]

    # Generate remaining states in sequence
    if (actual_seq_length > 1) {
      for (i in 2:actual_seq_length) {
        current_state_index <- match(sequence[i - 1], states_names)
        next_state_index <- sample(num_states_actual, 1, prob = trans_matrix[current_state_index, ])
        sequence[i] <- states_names[next_state_index]
      }
    }

    # Fill the row: first with valid states, then with NAs
    result_matrix[row, 1:actual_seq_length] <- sequence
    # The rest of the row will remain NA (initialized that way)
  }

  # Convert to data frame with standard column names
  result_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
  names(result_df) <- paste0("V", 1:seq_length)

  # Return based on include_params
  if (include_params) {
    return(list(
      sequences = result_df,
      transition_matrix = trans_matrix,
      initial_probabilities = init_probs,
      state_names = states_names
    ))
  }

  return(result_df)
}
