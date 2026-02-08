#' Simulate TNA Datasets (Sequences + Models + Probabilities)
#'
#' @description
#' Simulate complete TNA datasets including simulated sequences, fitted models,
#' and generating probabilities. This function combines the full simulation
#' workflow: probability generation -> sequence simulation -> model fitting.
#'
#' Use \code{\link{simulate_tna_networks}} if you only need the fitted models.
#'
#' Supports using realistic learning action verbs as state names for
#' educational research simulations.
#'
#' @param n_datasets Integer. Number of datasets to generate. Default: 1.
#' @param n_states Integer. Number of states in each network. Default: 5.
#' @param n_sequences Integer. Number of sequences to simulate per network.
#'   Default: 100.
#' @param seq_length Integer. Maximum length of each sequence. Default: 30.
#' @param states Character vector. Names for the states. If NULL and
#'   `use_learning_states = FALSE`, uses letters (A, B, C, ...).
#'   Ignored if `use_learning_states = TRUE`. Default: NULL.
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names (e.g., "Plan", "Monitor", "Read"). Default: TRUE.
#' @param categories Character vector. Which categories of learning
#'   verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", or "all".
#'   If NULL and `use_learning_states = TRUE`, a random category is selected.
#'   Only used if `use_learning_states = TRUE`. Default: NULL.
#' @param smart_select Logical. If TRUE, intelligently selects learning states
#'   based on n_states (small networks use fewer categories).
#'   Only used if `use_learning_states = TRUE`. Default: TRUE.
#' @param na_range Integer vector of length 2 (min, max) or single integer.
#'   Range of NA values per sequence. Default: c(0, 0).
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
#' @param include_data Logical. If TRUE, includes the generating
#'   sequence data in the output. Default: TRUE.
#' @param include_probs Logical. If TRUE, includes the generating
#'   transition matrix and initial probabilities in the output. Default: TRUE.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL,
#'   no seed is set. Default: NULL.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#'
#' @param state_names Deprecated. Use `states` instead.
#' @param learning_categories Deprecated. Use `categories` instead.
#' @param num_rows Deprecated. Use `n_sequences` instead.
#' @param max_seq_length Deprecated. Use `seq_length` instead.
#' @param min_na Deprecated. Use `na_range` instead.
#' @param max_na Deprecated. Use `na_range` instead.
#' @param include_probabilities Deprecated. Use `include_probs` instead.
#'
#' @return A list of length `n_datasets`. Each element is a list containing:
#' \describe{
#'   \item{model}{The fitted TNA model object.}
#'   \item{transition_probs}{The generating transition matrix
#'     (if `include_probs = TRUE`).}
#'   \item{initial_probs}{The generating initial probabilities
#'     (if `include_probs = TRUE`).}
#'   \item{sequences}{The simulated sequence data frame
#'     (if `include_data = TRUE`).}
#'   \item{params}{List of parameters used for this network.}
#' }
#'
#' @details
#' The function generates sequence data by:
#' \enumerate{
#'   \item Generating random transition and initial probabilities using
#'     `seqHMM::simulate_transition_probs()` and `seqHMM::simulate_initial_probs()`.
#'   \item Simulating sequences from those probabilities using either
#'     `simulate_sequences()` or `simulate_sequences_advanced()`.
#'   \item Fitting a TNA model to the sequences using `fit_network_model()`.
#' }
#'
#' **Learning States**: When `use_learning_states = TRUE`, state names are
#' drawn from a curated collection of 180+ student learning action verbs
#' organized into 7 categories:
#' \itemize{
#'   \item **metacognitive**: Plan, Monitor, Evaluate, Reflect, Regulate, ...
#'   \item **cognitive**: Read, Study, Analyze, Summarize, Memorize, ...
#'   \item **behavioral**: Practice, Annotate, Research, Review, Write, ...
#'   \item **social**: Collaborate, Discuss, Explain, Share, Teach, ...
#'   \item **motivational**: Focus, Persist, Explore, Strive, Commit, ...
#'   \item **affective**: Enjoy, Appreciate, Cope, Manage, Curious, ...
#'   \item **group_regulation**: Adapt, Cohesion, Consensus, Coregulate, ...
#' }
#'
#' @examples
#' \dontrun{
#' # Simplest usage - generates 1 dataset with sequences + model + probabilities
#' data <- simulate_tna_datasets(seed = 42)
#' data[[1]]$sequences        # The sequence data
#' data[[1]]$model            # The fitted TNA model
#' data[[1]]$transition_probs # The generating probabilities
#'
#' # Generate 5 datasets with learning states
#' data <- simulate_tna_datasets(
#'   n_datasets = 5,
#'   n_states = 4,
#'   seed = 42
#' )
#'
#' # Generate with specific learning category
#' learning_data <- simulate_tna_datasets(
#'   n_datasets = 5,
#'   n_states = 6,
#'   categories = c("metacognitive", "cognitive"),
#'   seed = 42
#' )
#'
#' # View the state names
#' learning_data[[1]]$params$state_names
#' # e.g., c("Plan", "Monitor", "Read", "Practice", "Discuss", "Focus")
#'
#' # Generate with letter names (disable learning states)
#' letter_data <- simulate_tna_datasets(
#'   n_datasets = 5,
#'   n_states = 8,
#'   use_learning_states = FALSE,
#'   seed = 123
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_tna_networks}} for generating TNA network objects,
#' \code{\link{get_learning_states}} for retrieving learning state verbs,
#' \code{\link{select_states}} for intelligent state selection,
#' \code{\link{list_learning_categories}} for viewing available categories,
#' \code{\link{simulate_sequences}} for basic sequence simulation,
#' \code{\link{fit_network_model}} for model fitting.
#'
#' @name simulate_tna_datasets
#' @rdname simulate_tna_datasets
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
simulate_tna_datasets <- function(n_datasets = 1,
                                   n_states = 8,
                                   n_sequences = 100,
                                   seq_length = 30,
                                   states = NULL,
                                   use_learning_states = TRUE,
                                   categories = NULL,
                                   smart_select = TRUE,
                                   na_range = c(0, 0),
                                   model_type = "tna",
                                   use_advanced = FALSE,
                                   stable_transitions = NULL,
                                   stability_prob = 0.95,
                                   unstable_mode = "random_jump",
                                   unstable_random_transition_prob = 0.5,
                                   include_data = TRUE,
                                   include_probs = TRUE,
                                   seed = NULL,
                                   verbose = TRUE,
                                   # Backward compatibility - old parameter names
                                   state_names = NULL,
                                   learning_categories = NULL,
                                   num_rows = NULL,
                                   max_seq_length = NULL,
                                   min_na = NULL,
                                   max_na = NULL,
                                   include_probabilities = NULL) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(state_names)) states <- state_names
  if (!is.null(learning_categories)) categories <- learning_categories
  if (!is.null(num_rows)) n_sequences <- num_rows
  if (!is.null(max_seq_length)) seq_length <- max_seq_length
  if (!is.null(include_probabilities)) include_probs <- include_probabilities

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

  # --- Input Validation ---
  stopifnot(
    is.numeric(n_datasets), length(n_datasets) == 1, n_datasets >= 1,
    is.numeric(n_states), length(n_states) == 1, n_states >= 2,
    is.numeric(n_sequences), length(n_sequences) == 1, n_sequences >= 1,
    is.numeric(seq_length), length(seq_length) == 1, seq_length >= 2,
    is.numeric(na_range), length(na_range) == 2, na_range[1] >= 0, na_range[2] >= na_range[1],
    is.character(model_type), length(model_type) == 1,
    model_type %in% c("tna", "ftna", "ctna", "atna"),
    is.logical(use_advanced), length(use_advanced) == 1,
    is.logical(use_learning_states), length(use_learning_states) == 1,
    is.logical(smart_select), length(smart_select) == 1,
    is.logical(include_data), length(include_data) == 1,
    is.logical(include_probs), length(include_probs) == 1,
    is.logical(verbose), length(verbose) == 1
  )

  # Set seed if provided (do this early for reproducible state selection)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Determine State Names ---
  if (use_learning_states) {
    # Random category selection if not specified
    if (is.null(categories)) {
      all_cats <- c("metacognitive", "cognitive", "behavioral",
                    "social", "motivational", "affective", "group_regulation")
      categories <- sample(all_cats, 1)
      if (verbose) {
        message(sprintf("Randomly selected learning category: %s", categories))
      }
    }

    if (smart_select && ("all" %in% categories || length(categories) > 2)) {
      # Use smart selection for balanced category representation
      states <- select_states(
        n_states = n_states,
        primary_categories = if ("all" %in% categories) NULL else categories
      )
    } else {
      # Direct selection from specified categories
      states <- get_learning_states(
        categories = categories,
        n = n_states
      )
    }
    if (verbose) {
      message(sprintf("Using learning states: %s", paste(states, collapse = ", ")))
    }
  } else if (is.null(states)) {
    # Default to letters or S1, S2, ...
    if (n_states <= 26) {
      states <- LETTERS[1:n_states]
    } else {
      states <- paste0("S", 1:n_states)
    }
  } else {
    stopifnot(is.character(states), length(states) >= n_states)
    states <- states[1:n_states]
  }

  # --- Generate Datasets ---
  datasets <- vector("list", n_datasets)

  for (i in seq_len(n_datasets)) {
    if (verbose) {
      message(sprintf("Generating dataset %d/%d...", i, n_datasets))
    }

    # Generate random transition and initial probabilities
    initial_probs <- seqHMM::simulate_initial_probs(n_states)
    names(initial_probs) <- states
    transition_probs <- seqHMM::simulate_transition_probs(n_states)
    dimnames(transition_probs) <- list(states, states)

    # Simulate sequences
    if (use_advanced) {
      sequences <- simulate_sequences_advanced(
        trans_matrix = transition_probs,
        init_probs = initial_probs,
        seq_length = seq_length,
        n_sequences = n_sequences,
        stable_transitions = stable_transitions,
        stability_prob = stability_prob,
        unstable_mode = unstable_mode,
        unstable_random_transition_prob = unstable_random_transition_prob,
        na_range = na_range,
        include_na = (na_range[2] > 0)
      )
    } else {
      sequences <- simulate_sequences(
        trans_matrix = transition_probs,
        init_probs = initial_probs,
        seq_length = seq_length,
        n_sequences = n_sequences,
        na_range = na_range,
        include_na = (na_range[2] > 0)
      )
    }

    # Fit TNA model
    model <- tryCatch(
      {
        fit_network_model(sequences, model_type)
      },
      error = function(e) {
        warning(sprintf("Failed to fit model for dataset %d: %s", i, e$message))
        NULL
      }
    )

    # Build result for this network
    result <- list(
      model = model,
      params = list(
        network_id = i,
        n_states = n_states,
        state_names = states,
        n_sequences = n_sequences,
        seq_length = seq_length,
        na_range = na_range,
        model_type = model_type,
        use_advanced = use_advanced,
        use_learning_states = use_learning_states,
        categories = if (use_learning_states) categories else NULL
      )
    )

    if (include_probs) {
      result$transition_probs <- transition_probs
      result$initial_probs <- initial_probs
    }

    if (include_data) {
      result$sequences <- sequences
    }

    datasets[[i]] <- result
  }

  # Name the list elements
  names(datasets) <- paste0("dataset_", seq_len(n_datasets))

  if (verbose) {
    successful <- sum(sapply(datasets, function(x) !is.null(x$model)))
    message(sprintf("Generated %d/%d datasets successfully.", successful, n_datasets))
  }

  return(datasets)
}


#' Simulate TNA Network Objects
#'
#' @description
#' Simulate multiple TNA network objects (fitted models) for simulation studies.
#' This function creates TNA model objects by simulating sequence data and
#' fitting models. For generating raw sequence data with associated probabilities,
#' use \code{\link{simulate_tna_datasets}}. For group TNA models, use
#' \code{\link{simulate_group_tna_networks}}.
#'
#' @param n_networks Integer. Number of networks to generate. Default: 1.
#' @param n_states Integer. Number of states/actions in each network. Default: 5.
#' @param n_sequences Integer. Number of sequences to simulate per network.
#'   Default: 100.
#' @param seq_length Integer. Maximum length of each sequence. Default: 30.
#' @param model_type Character. Type of TNA model to fit:
#'   "tna", "ftna", "ctna", "atna". Default: "tna".
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names. Default: TRUE.
#' @param categories Character vector or NULL. Which categories of
#'   learning verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", or "all".
#'   If NULL (default), randomly selects one category.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#' @param ... Additional parameters passed to \code{\link{simulate_tna_datasets}}.
#'
#' @param num_sequences Deprecated. Use `n_sequences` instead.
#' @param max_seq_length Deprecated. Use `seq_length` instead.
#' @param learning_categories Deprecated. Use `categories` instead.
#'
#' @return A list of length `n_networks` containing fitted TNA model objects.
#'   Each element is a tna model object (class "tna").
#'
#' @details
#' This function differs from \code{\link{simulate_tna_datasets}} in that it
#' returns only the fitted model objects, not the underlying sequence data or
#' generating probabilities. Use this when you need TNA network objects for
#' simulation studies or method comparisons.
#'
#' **Random Category Selection**: When `categories = NULL` and
#' `use_learning_states = TRUE`, a random learning category is selected.
#' Available categories: metacognitive, cognitive, behavioral, social,
#' motivational, affective, group_regulation.
#'
#' @examples
#' \dontrun{
#' # Generate 5 TNA networks with random learning category
#' nets <- simulate_tna_networks(n_networks = 5, seed = 42)
#'
#' # Generate networks with specific category
#' meta_nets <- simulate_tna_networks(
#'   n_networks = 3,
#'   n_states = 6,
#'   categories = "metacognitive",
#'   seed = 123
#' )
#'
#' # Generate filtered TNA networks
#' ftna_nets <- simulate_tna_networks(
#'   n_networks = 5,
#'   model_type = "ftna",
#'   seed = 456
#' )
#'
#' # Use the networks
#' plot(nets[[1]])
#' tna::centralities(nets[[1]])
#'
#' # Old parameter names still work
#' nets <- simulate_tna_networks(
#'   n_networks = 3,
#'   num_sequences = 50,
#'   max_seq_length = 25,
#'   seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_group_tna_networks}} for group TNA models,
#' \code{\link{simulate_tna_datasets}} for generating complete datasets with
#'   probabilities, \code{\link{fit_network_model}} for model fitting,
#' \code{\link{get_learning_states}} for learning state verbs.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
simulate_tna_networks <- function(n_networks = 1,
                                   n_states = 8,
                                   n_sequences = 100,
                                   seq_length = 30,
                                   model_type = "tna",
                                   use_learning_states = TRUE,
                                   categories = NULL,
                                   seed = NULL,
                                   verbose = TRUE,
                                   # Backward compatibility - old parameter names
                                   num_sequences = NULL,
                                   max_seq_length = NULL,
                                   learning_categories = NULL,
                                   ...) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(num_sequences)) n_sequences <- num_sequences
  if (!is.null(max_seq_length)) seq_length <- max_seq_length
  if (!is.null(learning_categories)) categories <- learning_categories

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Random category selection if not specified
  if (use_learning_states && is.null(categories)) {
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    categories <- sample(all_cats, 1)
    if (verbose) {
      message(sprintf("Randomly selected learning category: %s", categories))
    }
  }

  # Determine state names
  if (use_learning_states) {
    states <- get_learning_states(categories, n = n_states, seed = seed)
    if (verbose) {
      message(sprintf("Using learning states: %s", paste(states, collapse = ", ")))
    }
  } else {
    if (n_states <= 26) {
      states <- LETTERS[1:n_states]
    } else {
      states <- paste0("S", 1:n_states)
    }
  }

  networks <- list()

  for (i in seq_len(n_networks)) {
    if (verbose) {
      message(sprintf("Generating network %d/%d...", i, n_networks))
    }

    # Generate sequences and fit tna model
    seq_data <- simulate_tna_datasets(
      n_datasets = 1,
      n_states = n_states,
      states = states,
      use_learning_states = FALSE,
      n_sequences = n_sequences,
      seq_length = seq_length,
      model_type = model_type,
      seed = NULL,
      include_probs = FALSE,
      include_data = FALSE,
      verbose = FALSE,
      ...
    )
    networks[[paste0("network_", i)]] <- seq_data[[1]]$model
  }

  if (verbose) {
    successful <- sum(!sapply(networks, is.null))
    message(sprintf("Generated %d/%d networks successfully.", successful, n_networks))
  }

  return(networks)
}


#' Simulate Group TNA Network Objects
#'
#' @description
#' Simulate group TNA network objects (fitted group_tna models) for simulation
#' studies with grouped/clustered data. Each network contains multiple groups
#' (e.g., classrooms, teams) with their own transition patterns.
#'
#' @param n_groups Integer. Number of groups in the network. Default: 5.
#' @param n_actors Integer or integer vector. Number of actors per group.
#'   Accepts: single integer (fixed size), two integers like `c(8, 12)` (min/max),
#'   or a range like `5:15` (min/max taken from range). Default: 10.
#' @param n_states Integer. Number of states/actions in the network. Default: 5.
#' @param seq_length_range Integer vector of length 2. Range for sequence lengths
#'   per actor (min, max). Default: c(10, 30).
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names. Default: TRUE.
#' @param categories Character vector or NULL. Which categories of
#'   learning verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", or "all".
#'   If NULL (default), randomly selects one category.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#' @param ... Additional parameters passed to \code{\link{simulate_long_data}}.
#'
#' @param actors_per_group Deprecated. Use `n_actors` instead.
#' @param min_seq_length Deprecated. Use `seq_length_range` instead.
#' @param max_seq_length Deprecated. Use `seq_length_range` instead.
#' @param learning_categories Deprecated. Use `categories` instead.
#'
#' @return A group_tna model object (class "group_tna") containing:
#' \describe{
#'   \item{networks}{List of TNA networks, one per group}
#'   \item{data}{The underlying wide-format data}
#'   \item{group}{The grouping variable name}
#' }
#'
#' @details
#' This function generates a single group TNA network with multiple groups.
#' It simulates long-format sequence data with group structure, converts it
#' to wide format, and fits a group TNA model using \code{tna::group_model()}.
#'
#' **Use Cases**:
#' \itemize{
#'   \item Simulating classroom-level learning behavior data
#'   \item Generating team collaboration sequences
#'   \item Creating multi-group datasets for method comparison
#' }
#'
#' @examples
#' \dontrun{
#' # Generate a group TNA network with 4 groups, 15 actors each
#' group_net <- simulate_group_tna_networks(
#'   n_groups = 4,
#'   n_actors = 15,
#'   n_states = 5,
#'   seed = 42
#' )
#'
#' # Variable group sizes using range notation
#' var_net <- simulate_group_tna_networks(
#'   n_groups = 5,
#'   n_actors = c(8, 15),
#'   seq_length_range = c(5, 25),
#'   seed = 123
#' )
#'
#' # With specific learning category
#' ssrl_net <- simulate_group_tna_networks(
#'   n_groups = 6,
#'   n_actors = c(10, 20),
#'   categories = "group_regulation",
#'   seed = 456
#' )
#'
#' # Access individual group networks
#' names(group_net)
#' group_net[[1]]$weights
#'
#' # Old parameter names still work
#' group_net <- simulate_group_tna_networks(
#'   actors_per_group = 10,
#'   min_seq_length = 5,
#'   max_seq_length = 20,
#'   seed = 42
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_tna_networks}} for individual TNA models,
#' \code{\link{simulate_long_data}} for generating long-format group data,
#' \code{\link{fit_network_model}} for model fitting.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
simulate_group_tna_networks <- function(n_groups = 5,
                                         n_actors = 10,
                                         n_states = 8,
                                         seq_length_range = c(10, 30),
                                         use_learning_states = TRUE,
                                         categories = NULL,
                                         seed = NULL,
                                         verbose = TRUE,
                                         # Backward compatibility - old parameter names
                                         actors_per_group = NULL,
                                         min_seq_length = NULL,
                                         max_seq_length = NULL,
                                         learning_categories = NULL,
                                         ...) {
  # --- Backward compatibility: map old names to new names ---
  if (!is.null(actors_per_group)) n_actors <- actors_per_group
  if (!is.null(learning_categories)) categories <- learning_categories

  # Handle seq_length_range from min/max_seq_length
  if (!is.null(min_seq_length) || !is.null(max_seq_length)) {
    min_len <- if (!is.null(min_seq_length)) min_seq_length else seq_length_range[1]
    max_len <- if (!is.null(max_seq_length)) max_seq_length else seq_length_range[2]
    seq_length_range <- c(min_len, max_len)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Handle n_actors: support single value, c(min, max), or min:max
  if (length(n_actors) == 1) {
    actors_range <- c(n_actors, n_actors)
  } else {
    actors_range <- c(min(n_actors), max(n_actors))
  }

  # Random category selection if not specified
  if (use_learning_states && is.null(categories)) {
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    categories <- sample(all_cats, 1)
    if (verbose) {
      message(sprintf("Randomly selected learning category: %s", categories))
    }
  }

  # Determine state names
  if (use_learning_states) {
    states <- get_learning_states(categories, n = n_states, seed = seed)
    if (verbose) {
      message(sprintf("Using learning states: %s", paste(states, collapse = ", ")))
    }
  } else {
    if (n_states <= 26) {
      states <- LETTERS[1:n_states]
    } else {
      states <- paste0("S", 1:n_states)
    }
  }

  if (verbose) {
    if (actors_range[1] == actors_range[2]) {
      actors_msg <- sprintf("%d actors/group", actors_range[1])
    } else {
      actors_msg <- sprintf("%d-%d actors/group", actors_range[1], actors_range[2])
    }
    message(sprintf("Generating group TNA network: %d groups, %s, seq length %d-%d",
                    n_groups, actors_msg, seq_length_range[1], seq_length_range[2]))
  }

  # Generate grouped long-format data
  long_data <- simulate_long_data(
    n_groups = n_groups,
    n_actors = actors_range,
    states = states,
    seq_length_range = seq_length_range,
    seed = NULL,
    ...
  )

  # Convert to wide format
  wide_data <- long_to_wide(long_data, id_col = "Actor")

  # Add Group column back to wide_data for group_model
  actor_group_map <- unique(long_data[, c("Actor", "Group")])
  wide_data <- merge(wide_data, actor_group_map, by = "Actor", all.x = TRUE)

  # Remove Actor column - group_model only needs sequence columns + group
  wide_data$Actor <- NULL

  # Fit group model
  model <- tryCatch(
    {
      fit_network_model(wide_data, "group_tna", group = "Group")
    },
    error = function(e) {
      stop(sprintf("Failed to fit group model: %s", e$message))
    }
  )

  if (verbose) {
    message("Group TNA network generated successfully.")
  }

  return(model)
}


#' Simulate TNA Transition Matrix with Node Groupings
#'
#' @description
#' Simulate a transition matrix with node groupings, compatible with
#' `tna::plot_htna()` (hierarchical) and `tna::plot_mlna()` (multilevel)
#' visualizations. By default creates a 25-node matrix (5 nodes x 5 types)
#' using learning category names.
#'
#' This is a convenience wrapper around \code{\link{simulate_htna}}.
#'
#' @param nodes_per_group Integer. Number of nodes per group. Default: 5.
#' @param group_names Character vector. Names for each group. Default uses
#'   learning categories: "Metacognitive", "Cognitive", "Behavioral",
#'   "Social", "Motivational".
#' @param n_groups Integer. Number of groups. Default: 5.
#' @param edge_prob_range Numeric vector of length 2. Range for edge weights
#'   `c(min, max)`. Default: c(0, 1).
#' @param self_loops Logical. Allow self-loops (diagonal elements). Default: FALSE.
#' @param use_learning_states Logical. Use learning state verbs as node names.
#'   Default: TRUE.
#' @param categories Character vector. Categories for node names, one per group.
#'   Default: c("metacognitive", "cognitive", "behavioral", "social", "motivational").
#' @param within_prob Numeric. Probability of edges within each group. Default: 0.4.
#' @param between_prob Numeric. Probability of edges between groups. Default: 0.15.
#' @param node_prefix Character. Prefix for node names when not using learning
#'   states. Default: "N".
#' @param seed Integer or NULL. Random seed. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @param learning_categories Deprecated. Use `categories` instead.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{matrix}{Square transition matrix (rows sum to 1) with named rows/columns.}
#'   \item{node_types}{Named list mapping group names to node names.
#'     Use as `node_types` for `plot_htna()` or as `layers` for `plot_mlna()`.}
#' }
#'
#' @details
#' This function generates a random transition matrix and node groupings.
#' The output can be used with both hierarchical and multilevel TNA plots:
#'
#' \itemize{
#'   \item For `plot_htna()`: use `net$node_types` directly
#'   \item For `plot_mlna()`: use `layers = net$node_types`
#' }
#'
#' @examples
#' \dontrun{
#' # Default: 5 groups x 5 nodes = 25 node matrix
#' net <- simulate_tna_matrix(seed = 42)
#' net$matrix
#' net$node_types  # Metacognitive, Cognitive, Behavioral, Social, Motivational
#'
#' # Use with plot_htna
#' plot_htna(net$matrix, net$node_types, layout = "polygon")
#'
#' # Use with plot_mlna
#' plot_mlna(net$matrix, layers = net$node_types)
#'
#' # Custom group names (3 groups)
#' net <- simulate_tna_matrix(
#'   nodes_per_group = 6,
#'   group_names = c("Macro", "Meso", "Micro"),
#'   seed = 42
#' )
#'
#' # Custom categories per group
#' net <- simulate_tna_matrix(
#'   nodes_per_group = 4,
#'   group_names = c("Teacher", "Student", "System"),
#'   categories = c("metacognitive", "cognitive", "behavioral"),
#'   seed = 123
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_htna}} for the underlying function,
#' \code{\link{simulate_matrix}} for basic matrix simulation,
#' \code{\link{simulate_tna_networks}} for TNA model objects.
#'
#' @export
simulate_tna_matrix <- function(nodes_per_group = 5,
                                 group_names = c("Metacognitive", "Cognitive",
                                                 "Behavioral", "Social", "Motivational"),
                                 n_groups = 5,
                                 edge_prob_range = c(0, 1),
                                 self_loops = FALSE,
                                 use_learning_states = TRUE,
                                 categories = c("metacognitive", "cognitive",
                                                "behavioral", "social", "motivational"),
                                 within_prob = 0.4,
                                 between_prob = 0.15,
                                 node_prefix = "N",
                                 seed = NULL,
                                 verbose = TRUE,
                                 # Backward compatibility
                                 learning_categories = NULL) {
  # Backward compatibility
  if (!is.null(learning_categories)) categories <- learning_categories

  # Determine number of groups from group_names if custom names provided
  if (!identical(group_names, c("Metacognitive", "Cognitive",
                                "Behavioral", "Social", "Motivational"))) {
    n_groups <- length(group_names)
  }

  # Adjust categories to match n_groups
  if (length(categories) != n_groups) {
    default_cats <- c("metacognitive", "cognitive", "behavioral",
                      "social", "motivational", "affective", "group_regulation")
    if (n_groups <= length(default_cats)) {
      categories <- default_cats[1:n_groups]
    } else {
      categories <- c(default_cats, rep("all", n_groups - length(default_cats)))
    }
  }

  # Adjust group_names if n_groups changed but group_names wasn't custom
  if (n_groups != 5 && identical(group_names, c("Metacognitive", "Cognitive",
                                                 "Behavioral", "Social", "Motivational"))) {
    default_names <- c("Metacognitive", "Cognitive", "Behavioral",
                       "Social", "Motivational", "Affective", "GroupRegulation")
    if (n_groups <= length(default_names)) {
      group_names <- default_names[1:n_groups]
    } else {
      group_names <- c(default_names, paste0("Type", (length(default_names) + 1):n_groups))
    }
  }

  # Call simulate_htna
  result <- simulate_htna(
    n_nodes = nodes_per_group,
    n_types = n_groups,
    type_names = group_names,
    within_prob = within_prob,
    between_prob = between_prob,
    weight_range = edge_prob_range,
    allow_self_loops = self_loops,
    categories = categories,
    seed = seed
  )

  if (verbose) {
    total_nodes <- nodes_per_group * n_groups
    message(sprintf("Generated TNA matrix: %d groups, %d total nodes",
                    n_groups, total_nodes))
    for (i in seq_len(n_groups)) {
      message(sprintf("  %s: %s", group_names[i],
                      paste(result$node_types[[group_names[i]]], collapse = ", ")))
    }
  }

  return(list(
    matrix = result$matrix,
    node_types = result$node_types
  ))
}
