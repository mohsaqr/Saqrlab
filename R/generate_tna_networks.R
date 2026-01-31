#' Generate TNA Networks from Simulated Sequences
#'
#' @description
#' Generate multiple TNA network objects by simulating Markov chain sequences
#' and fitting TNA models. This function combines the simulation workflow:
#' probability generation -> sequence simulation -> model fitting.
#'
#' Supports using realistic learning action verbs as state names for
#' educational research simulations.
#'
#' @param n_networks Integer. Number of networks to generate. Default: 10.
#' @param n_states Integer. Number of states in each network. Default: 4.
#' @param state_names Character vector. Names for the states. If NULL and
#'   `use_learning_states = FALSE`, uses letters (A, B, C, ...).
#'   Ignored if `use_learning_states = TRUE`. Default: NULL.
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names (e.g., "Plan", "Monitor", "Read"). Default: FALSE.
#' @param learning_categories Character vector. Which categories of learning
#'   verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", or "all".
#'   Only used if `use_learning_states = TRUE`. Default: "all".
#' @param smart_select Logical. If TRUE, intelligently selects learning states
#'   based on n_states (small networks use fewer categories).
#'   Only used if `use_learning_states = TRUE`. Default: TRUE.
#' @param num_rows Integer. Number of sequences to simulate per network.
#'   Default: 100.
#' @param max_seq_length Integer. Maximum length of each sequence. Default: 30.
#' @param min_na Integer. Minimum NAs per sequence. Default: 0.
#' @param max_na Integer. Maximum NAs per sequence. Default: 0.
#' @param model_type Character. Type of TNA model to fit. One of:
#'   "tna", "ftna", "ctna", "atna". Default: "tna".
#' @param use_advanced Logical. If TRUE, uses `simulate_sequences_advanced()`
#'   with stability modes. If FALSE, uses basic `simulate_sequences()`.
#'   Default: FALSE.
#' @param stable_transitions List of character vectors defining stable
#'   transitions (only used if `use_advanced = TRUE`). Default: NULL.
#' @param stability_prob Numeric (0 to 1). Probability of following stable
#'   transitions (only used if `use_advanced = TRUE`). Default: 0.95.
#' @param unstable_mode Character. Mode for unstable transitions:
#'   "random_jump", "perturb_prob", or "unlikely_jump".
#'   (Only used if `use_advanced = TRUE`). Default: "random_jump".
#' @param unstable_random_transition_prob Numeric (0 to 1). Probability of
#'   unstable action (only used if `use_advanced = TRUE`). Default: 0.5.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL,
#'   no seed is set. Default: NULL.
#' @param include_probabilities Logical. If TRUE, includes the generating
#'   transition matrix and initial probabilities in the output. Default: TRUE.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#'
#' @return A list of length `n_networks`. Each element is a list containing:
#' \describe{
#'   \item{model}{The fitted TNA model object.}
#'   \item{transition_probs}{The generating transition matrix
#'     (if `include_probabilities = TRUE`).}
#'   \item{initial_probs}{The generating initial probabilities
#'     (if `include_probabilities = TRUE`).}
#'   \item{sequences}{The simulated sequence data frame
#'     (if `include_probabilities = TRUE`).}
#'   \item{params}{List of parameters used for this network.}
#' }
#'
#' @details
#' The function generates networks by:
#' \enumerate{
#'   \item Generating random transition and initial probabilities using
#'     `seqHMM::simulate_transition_probs()` and `seqHMM::simulate_initial_probs()`.
#'   \item Simulating sequences from those probabilities using either
#'     `simulate_sequences()` or `simulate_sequences_advanced()`.
#'   \item Fitting a TNA model to the sequences using `fit_network_model()`.
#' }
#'
#' **Learning States**: When `use_learning_states = TRUE`, state names are
#' drawn from a curated collection of 180 student learning action verbs
#' organized into 6 categories:
#' \itemize{
#'   \item **metacognitive**: Plan, Monitor, Evaluate, Reflect, Regulate, ...
#'   \item **cognitive**: Read, Study, Analyze, Summarize, Memorize, ...
#'   \item **behavioral**: Practice, Annotate, Research, Review, Write, ...
#'   \item **social**: Collaborate, Discuss, Explain, Share, Teach, ...
#'   \item **motivational**: Focus, Persist, Explore, Strive, Commit, ...
#'   \item **affective**: Enjoy, Appreciate, Cope, Manage, Curious, ...
#' }
#'
#' @examples
#' \dontrun{
#' # Generate 5 simple TNA networks with letter names
#' networks <- generate_tna_networks(
#'   n_networks = 5,
#'   n_states = 4,
#'   seed = 42
#' )
#'
#' # Generate networks with learning action verbs
#' learning_nets <- generate_tna_networks(
#'   n_networks = 5,
#'   n_states = 6,
#'   use_learning_states = TRUE,
#'   seed = 42
#' )
#'
#' # View the state names
#' learning_nets[[1]]$params$state_names
#' # e.g., c("Plan", "Monitor", "Read", "Practice", "Discuss", "Focus")
#'
#' # Focus on metacognitive and cognitive states only
#' srl_networks <- generate_tna_networks(
#'   n_networks = 5,
#'   n_states = 8,
#'   use_learning_states = TRUE,
#'   learning_categories = c("metacognitive", "cognitive"),
#'   seed = 123
#' )
#'
#' # Social learning focus
#' social_nets <- generate_tna_networks(
#'   n_networks = 3,
#'   n_states = 6,
#'   use_learning_states = TRUE,
#'   learning_categories = c("social", "behavioral"),
#'   seed = 456
#' )
#'
#' # Large network with all categories
#' big_network <- generate_tna_networks(
#'   n_networks = 1,
#'   n_states = 20,
#'   use_learning_states = TRUE,
#'   learning_categories = "all",
#'   seed = 789
#' )
#'
#' # Advanced simulation with stability modes
#' networks_adv <- generate_tna_networks(
#'   n_networks = 3,
#'   n_states = 5,
#'   use_learning_states = TRUE,
#'   learning_categories = "metacognitive",
#'   use_advanced = TRUE,
#'   stability_prob = 0.9,
#'   unstable_mode = "unlikely_jump",
#'   seed = 123
#' )
#' }
#'
#' @seealso
#' \code{\link{get_learning_states}} for retrieving learning state verbs,
#' \code{\link{smart_select_states}} for intelligent state selection,
#' \code{\link{list_learning_categories}} for viewing available categories,
#' \code{\link{simulate_sequences}} for basic sequence simulation,
#' \code{\link{fit_network_model}} for model fitting.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
generate_tna_networks <- function(n_networks = 10,
                                   n_states = 4,
                                   state_names = NULL,
                                   use_learning_states = FALSE,
                                   learning_categories = "all",
                                   smart_select = TRUE,
                                   num_rows = 100,
                                   max_seq_length = 30,
                                   min_na = 0,
                                   max_na = 0,
                                   model_type = "tna",
                                   use_advanced = FALSE,
                                   stable_transitions = NULL,
                                   stability_prob = 0.95,
                                   unstable_mode = "random_jump",
                                   unstable_random_transition_prob = 0.5,
                                   seed = NULL,
                                   include_probabilities = TRUE,
                                   verbose = TRUE) {
  # --- Input Validation ---
  stopifnot(
    is.numeric(n_networks), length(n_networks) == 1, n_networks >= 1,
    is.numeric(n_states), length(n_states) == 1, n_states >= 2,
    is.numeric(num_rows), length(num_rows) == 1, num_rows >= 1,
    is.numeric(max_seq_length), length(max_seq_length) == 1, max_seq_length >= 2,
    is.numeric(min_na), length(min_na) == 1, min_na >= 0,
    is.numeric(max_na), length(max_na) == 1, max_na >= min_na,
    is.character(model_type), length(model_type) == 1,
    model_type %in% c("tna", "ftna", "ctna", "atna"),
    is.logical(use_advanced), length(use_advanced) == 1,
    is.logical(use_learning_states), length(use_learning_states) == 1,
    is.logical(smart_select), length(smart_select) == 1,
    is.logical(include_probabilities), length(include_probabilities) == 1,
    is.logical(verbose), length(verbose) == 1
  )

  # Set seed if provided (do this early for reproducible state selection)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Determine State Names ---
  if (use_learning_states) {
    if (smart_select && ("all" %in% learning_categories || length(learning_categories) > 2)) {
      # Use smart selection for balanced category representation
      state_names <- smart_select_states(
        n_states = n_states,
        primary_categories = if ("all" %in% learning_categories) NULL else learning_categories
      )
    } else {
      # Direct selection from specified categories
      state_names <- get_learning_states(
        categories = learning_categories,
        n = n_states
      )
    }
    if (verbose) {
      message(sprintf("Using learning states: %s", paste(state_names, collapse = ", ")))
    }
  } else if (is.null(state_names)) {
    # Default to letters or S1, S2, ...
    if (n_states <= 26) {
      state_names <- LETTERS[1:n_states]
    } else {
      state_names <- paste0("S", 1:n_states)
    }
  } else {
    stopifnot(is.character(state_names), length(state_names) >= n_states)
    state_names <- state_names[1:n_states]
  }

  # --- Generate Networks ---
  networks <- vector("list", n_networks)

  for (i in seq_len(n_networks)) {
    if (verbose) {
      message(sprintf("Generating network %d/%d...", i, n_networks))
    }

    # Generate random transition and initial probabilities
    initial_probs <- seqHMM::simulate_initial_probs(n_states)
    names(initial_probs) <- state_names
    transition_probs <- seqHMM::simulate_transition_probs(n_states)
    dimnames(transition_probs) <- list(state_names, state_names)

    # Simulate sequences
    if (use_advanced) {
      sequences <- simulate_sequences_advanced(
        transition_matrix = transition_probs,
        initial_probabilities = initial_probs,
        max_seq_length = max_seq_length,
        num_rows = num_rows,
        stable_transitions = stable_transitions,
        stability_prob = stability_prob,
        unstable_mode = unstable_mode,
        unstable_random_transition_prob = unstable_random_transition_prob,
        min_na = min_na,
        max_na = max_na,
        include_na = (max_na > 0)
      )
    } else {
      sequences <- simulate_sequences(
        transition_matrix = transition_probs,
        initial_probabilities = initial_probs,
        max_seq_length = max_seq_length,
        num_rows = num_rows,
        min_na = min_na,
        max_na = max_na,
        include_na = (max_na > 0)
      )
    }

    # Fit TNA model
    model <- tryCatch(
      {
        fit_network_model(sequences, model_type)
      },
      error = function(e) {
        warning(sprintf("Failed to fit model for network %d: %s", i, e$message))
        NULL
      }
    )

    # Build result for this network
    result <- list(
      model = model,
      params = list(
        network_id = i,
        n_states = n_states,
        state_names = state_names,
        num_rows = num_rows,
        max_seq_length = max_seq_length,
        min_na = min_na,
        max_na = max_na,
        model_type = model_type,
        use_advanced = use_advanced,
        use_learning_states = use_learning_states,
        learning_categories = if (use_learning_states) learning_categories else NULL
      )
    )

    if (include_probabilities) {
      result$transition_probs <- transition_probs
      result$initial_probs <- initial_probs
      result$sequences <- sequences
    }

    networks[[i]] <- result
  }

  # Name the list elements
  names(networks) <- paste0("network_", seq_len(n_networks))

  if (verbose) {
    successful <- sum(sapply(networks, function(x) !is.null(x$model)))
    message(sprintf("Generated %d/%d networks successfully.", successful, n_networks))
  }

  return(networks)
}
