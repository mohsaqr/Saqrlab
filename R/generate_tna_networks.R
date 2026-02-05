#' Generate Sequence Data with Fitted TNA Models
#'
#' @description
#' Generate multiple simulated sequence datasets by simulating Markov chain
#' sequences and fitting TNA models. This function combines the simulation
#' workflow: probability generation -> sequence simulation -> model fitting.
#'
#' Supports using realistic learning action verbs as state names for
#' educational research simulations.
#'
#' @param n_networks Integer. Number of sequence datasets to generate. Default: 10.
#' @param n_states Integer. Number of states in each network. Default: 4.
#' @param state_names Character vector. Names for the states. If NULL and
#'   `use_learning_states = FALSE`, uses letters (A, B, C, ...).
#'   Ignored if `use_learning_states = TRUE`. Default: NULL.
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names (e.g., "Plan", "Monitor", "Read"). Default: FALSE.
#' @param learning_categories Character vector. Which categories of learning
#'   verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", or "all".
#'   If NULL and `use_learning_states = TRUE`, a random category is selected.
#'   Only used if `use_learning_states = TRUE`. Default: NULL.
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
#' # Generate 5 simple sequence datasets with letter names
#' data <- generate_sequence_data(
#'   n_networks = 5,
#'   n_states = 4,
#'   seed = 42
#' )
#'
#' # Generate data with learning action verbs
#' learning_data <- generate_sequence_data(
#'   n_networks = 5,
#'   n_states = 6,
#'   use_learning_states = TRUE,
#'   seed = 42
#' )
#'
#' # View the state names
#' learning_data[[1]]$params$state_names
#' # e.g., c("Plan", "Monitor", "Read", "Practice", "Discuss", "Focus")
#'
#' # Focus on metacognitive and cognitive states only
#' srl_data <- generate_sequence_data(
#'   n_networks = 5,
#'   n_states = 8,
#'   use_learning_states = TRUE,
#'   learning_categories = c("metacognitive", "cognitive"),
#'   seed = 123
#' )
#'
#' # Random category selection (default when learning_categories = NULL)
#' random_cat_data <- generate_sequence_data(
#'   n_networks = 3,
#'   n_states = 5,
#'   use_learning_states = TRUE,
#'   seed = 456
#' )
#'
#' # Large dataset with all categories
#' big_data <- generate_sequence_data(
#'   n_networks = 1,
#'   n_states = 20,
#'   use_learning_states = TRUE,
#'   learning_categories = "all",
#'   seed = 789
#' )
#'
#' # Advanced simulation with stability modes
#' data_adv <- generate_sequence_data(
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
#' \code{\link{generate_tna_networks}} for generating TNA network objects,
#' \code{\link{get_learning_states}} for retrieving learning state verbs,
#' \code{\link{smart_select_states}} for intelligent state selection,
#' \code{\link{list_learning_categories}} for viewing available categories,
#' \code{\link{simulate_sequences}} for basic sequence simulation,
#' \code{\link{fit_network_model}} for model fitting.
#'
#' @name generate_sequence_data
#' @rdname generate_sequence_data
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
generate_sequence_data <- function(n_networks = 10,
                                   n_states = 4,
                                   state_names = NULL,
                                   use_learning_states = FALSE,
                                   learning_categories = NULL,
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
    # Random category selection if not specified
    if (is.null(learning_categories)) {
      all_cats <- c("metacognitive", "cognitive", "behavioral",
                    "social", "motivational", "affective", "group_regulation")
      learning_categories <- sample(all_cats, 1)
      if (verbose) {
        message(sprintf("Randomly selected learning category: %s", learning_categories))
      }
    }

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


#' @rdname generate_sequence_data
#' @export
generate_tna_data <- generate_sequence_data


#' Generate TNA Network Objects
#'
#' @description
#' Generate multiple TNA network objects (fitted models) for simulation studies.
#' This function creates TNA model objects by simulating sequence data and
#' fitting models. For generating raw sequence data with associated probabilities,
#' use \code{\link{generate_sequence_data}}. For group TNA models, use
#' \code{\link{generate_group_tna_networks}}.
#'
#' @param n_networks Integer. Number of networks to generate. Default: 10.
#' @param model_type Character. Type of TNA model to fit:
#'   "tna", "ftna", "ctna", "atna". Default: "tna".
#' @param n_states Integer. Number of states/actions in each network. Default: 5.
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names. Default: TRUE.
#' @param learning_categories Character vector or NULL. Which categories of
#'   learning verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", or "all".
#'   If NULL (default), randomly selects one category.
#' @param num_sequences Integer. Number of sequences to simulate per network.
#'   Default: 100.
#' @param max_seq_length Integer. Maximum length of each sequence. Default: 30.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#' @param ... Additional parameters passed to \code{\link{generate_sequence_data}}.
#'
#' @return A list of length `n_networks` containing fitted TNA model objects.
#'   Each element is a tna model object (class "tna").
#'
#' @details
#' This function differs from \code{\link{generate_sequence_data}} in that it
#' returns only the fitted model objects, not the underlying sequence data or
#' generating probabilities. Use this when you need TNA network objects for
#' simulation studies or method comparisons.
#'
#' **Random Category Selection**: When `learning_categories = NULL` and
#' `use_learning_states = TRUE`, a random learning category is selected.
#' Available categories: metacognitive, cognitive, behavioral, social,
#' motivational, affective, group_regulation.
#'
#' @examples
#' \dontrun{
#' # Generate 5 TNA networks with random learning category
#' nets <- generate_tna_networks(n_networks = 5, seed = 42)
#'
#' # Generate networks with specific category
#' meta_nets <- generate_tna_networks(
#'   n_networks = 3,
#'   n_states = 6,
#'   learning_categories = "metacognitive",
#'   seed = 123
#' )
#'
#' # Generate filtered TNA networks
#' ftna_nets <- generate_tna_networks(
#'   n_networks = 5,
#'   model_type = "ftna",
#'   seed = 456
#' )
#'
#' # Use the networks
#' plot(nets[[1]])
#' tna::centralities(nets[[1]])
#' }
#'
#' @seealso
#' \code{\link{generate_group_tna_networks}} for group TNA models,
#' \code{\link{generate_sequence_data}} for generating sequence data with
#'   probabilities, \code{\link{fit_network_model}} for model fitting,
#' \code{\link{get_learning_states}} for learning state verbs.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
generate_tna_networks <- function(n_networks = 10,
                                   model_type = "tna",
                                   n_states = 5,
                                   use_learning_states = TRUE,
                                   learning_categories = NULL,
                                   num_sequences = 100,
                                   max_seq_length = 30,
                                   seed = NULL,
                                   verbose = TRUE,
                                   ...) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Random category selection if not specified
  if (use_learning_states && is.null(learning_categories)) {
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    learning_categories <- sample(all_cats, 1)
    if (verbose) {
      message(sprintf("Randomly selected learning category: %s", learning_categories))
    }
  }

  # Determine state names
  if (use_learning_states) {
    state_names <- get_learning_states(learning_categories, n = n_states, seed = seed)
    if (verbose) {
      message(sprintf("Using learning states: %s", paste(state_names, collapse = ", ")))
    }
  } else {
    if (n_states <= 26) {
      state_names <- LETTERS[1:n_states]
    } else {
      state_names <- paste0("S", 1:n_states)
    }
  }

  networks <- list()

  for (i in seq_len(n_networks)) {
    if (verbose) {
      message(sprintf("Generating network %d/%d...", i, n_networks))
    }

    # Generate sequences and fit tna model
    seq_data <- generate_sequence_data(
      n_networks = 1,
      n_states = n_states,
      state_names = state_names,
      use_learning_states = FALSE,
      num_rows = num_sequences,
      max_seq_length = max_seq_length,
      model_type = model_type,
      seed = NULL,
      include_probabilities = FALSE,
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


#' Generate Group TNA Network Objects
#'
#' @description
#' Generate group TNA network objects (fitted group_tna models) for simulation
#' studies with grouped/clustered data. Each network contains multiple groups
#' (e.g., classrooms, teams) with their own transition patterns.
#'
#' @param n_groups Integer. Number of groups in the network. Default: 5.
#' @param actors_per_group Integer or integer vector. Number of actors per group.
#'   Accepts: single integer (fixed size), two integers like `c(8, 12)` (min/max),
#'   or a range like `5:15` (min/max taken from range). Default: 10.
#' @param n_states Integer. Number of states/actions in the network. Default: 5.
#' @param min_seq_length Integer. Minimum sequence length per actor. Default: 10.
#' @param max_seq_length Integer. Maximum sequence length per actor. Default: 30.
#' @param use_learning_states Logical. If TRUE, uses learning action verbs
#'   as state names. Default: TRUE.
#' @param learning_categories Character vector or NULL. Which categories of
#'   learning verbs to use. Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", or "all".
#'   If NULL (default), randomly selects one category.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#' @param ... Additional parameters passed to \code{\link{simulate_long_data}}.
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
#' group_net <- generate_group_tna_networks(
#'   n_groups = 4,
#'   actors_per_group = 15,
#'   n_states = 5,
#'   seed = 42
#' )
#'
#' # Variable group sizes using range notation
#' var_net <- generate_group_tna_networks(
#'   n_groups = 5,
#'   actors_per_group = 8:15,
#'   min_seq_length = 5,
#'   max_seq_length = 25,
#'   seed = 123
#' )
#'
#' # With specific learning category
#' ssrl_net <- generate_group_tna_networks(
#'   n_groups = 6,
#'   actors_per_group = c(10, 20),
#'   learning_categories = "group_regulation",
#'   seed = 456
#' )
#'
#' # Access individual group networks
#' names(group_net)
#' group_net[[1]]$weights
#' }
#'
#' @seealso
#' \code{\link{generate_tna_networks}} for individual TNA models,
#' \code{\link{simulate_long_data}} for generating long-format group data,
#' \code{\link{fit_network_model}} for model fitting.
#'
#' @importFrom seqHMM simulate_initial_probs simulate_transition_probs
#' @import tna
#' @export
generate_group_tna_networks <- function(n_groups = 5,
                                         actors_per_group = 10,
                                         n_states = 5,
                                         min_seq_length = 10,
                                         max_seq_length = 30,
                                         use_learning_states = TRUE,
                                         learning_categories = NULL,
                                         seed = NULL,
                                         verbose = TRUE,
                                         ...) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Handle actors_per_group: support single value, c(min, max), or min:max

  if (length(actors_per_group) == 1) {
    actors_range <- c(actors_per_group, actors_per_group)
  } else {
    actors_range <- c(min(actors_per_group), max(actors_per_group))
  }

  # Random category selection if not specified
  if (use_learning_states && is.null(learning_categories)) {
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    learning_categories <- sample(all_cats, 1)
    if (verbose) {
      message(sprintf("Randomly selected learning category: %s", learning_categories))
    }
  }

  # Determine state names
  if (use_learning_states) {
    state_names <- get_learning_states(learning_categories, n = n_states, seed = seed)
    if (verbose) {
      message(sprintf("Using learning states: %s", paste(state_names, collapse = ", ")))
    }
  } else {
    if (n_states <= 26) {
      state_names <- LETTERS[1:n_states]
    } else {
      state_names <- paste0("S", 1:n_states)
    }
  }

  if (verbose) {
    if (actors_range[1] == actors_range[2]) {
      actors_msg <- sprintf("%d actors/group", actors_range[1])
    } else {
      actors_msg <- sprintf("%d-%d actors/group", actors_range[1], actors_range[2])
    }
    message(sprintf("Generating group TNA network: %d groups, %s, seq length %d-%d",
                    n_groups, actors_msg, min_seq_length, max_seq_length))
  }

  # Generate grouped long-format data
  long_data <- simulate_long_data(
    n_groups = n_groups,
    actors_per_group = actors_range,
    actions = state_names,
    seq_length_range = c(min_seq_length, max_seq_length),
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


#' Generate TNA Transition Matrix with Node Groupings
#'
#' @description
#' Generate a transition matrix with node groupings, compatible with
#' `tna::plot_htna()` (hierarchical) and `tna::plot_mlna()` (multilevel)
#' visualizations.
#'
#' @param nodes_per_group Integer or integer vector. Number of nodes per group.
#'   Accepts: single integer (same for all), vector matching n_groups,
#'   or range like `4:7` (random per group). Default: 5.
#' @param group_names Character vector. Names for each group (e.g.,
#'   `c("Teacher", "Student", "System")` or `c("Macro", "Meso", "Micro")`).
#'   If NULL, generates default names. Default: NULL.
#' @param n_groups Integer. Number of groups. Only used if `group_names` is NULL
#'   and `nodes_per_group` is a single integer. Default: 3.
#' @param edge_prob_range Numeric vector of length 2. Range for edge weights
#'   `c(min, max)`. Default: c(0, 0.3).
#' @param self_loops Logical. Allow self-loops (diagonal elements). Default: FALSE.
#' @param use_learning_states Logical. Use learning state verbs as node names.
#'   Default: TRUE.
#' @param learning_categories Character vector or NULL. Categories for node names.
#'   If NULL with `use_learning_states = TRUE`, uses one random category per group.
#'   Can be single (same for all) or one per group. Default: NULL.
#' @param node_prefix Character. Prefix for node names when not using learning
#'   states. Default: "N".
#' @param seed Integer or NULL. Random seed. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{matrix}{Square transition matrix with named rows/columns}
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
#' # Basic: 3 groups, 5 nodes each
#' net <- generate_tna_matrix(nodes_per_group = 5, n_groups = 3, seed = 42)
#'
#' # Use with plot_htna
#' plot_htna(net$matrix, net$node_types, layout = "polygon")
#' plot_htna(net$matrix, net$node_types, layout = "circular")
#'
#' # Use with plot_mlna (layers = node_types)
#' plot_mlna(net$matrix, layers = net$node_types, layout = "spring")
#'
#' # Custom group names
#' net <- generate_tna_matrix(
#'   nodes_per_group = c(5, 5, 5),
#'   group_names = c("Teacher", "Student", "System"),
#'   seed = 42
#' )
#'
#' # For multilevel with layer names
#' net <- generate_tna_matrix(
#'   nodes_per_group = 7,
#'   group_names = c("Macro", "Meso", "Micro"),
#'   learning_categories = c("metacognitive", "cognitive", "behavioral"),
#'   seed = 123
#' )
#' plot_mlna(net$matrix, layers = net$node_types, minimum = 0.18)
#'
#' # 4 groups for rectangle layout
#' net <- generate_tna_matrix(
#'   nodes_per_group = c(3, 2, 3, 3),
#'   group_names = c("Input", "Process", "Output", "Storage"),
#'   seed = 456
#' )
#' plot_htna(net$matrix, net$node_types)  # Auto-detects rectangle
#'
#' # Variable group sizes with range
#' net <- generate_tna_matrix(nodes_per_group = 4:8, n_groups = 3, seed = 789)
#' }
#'
#' @seealso
#' \code{\link{generate_tna_networks}} for TNA model objects,
#' \code{\link{generate_group_tna_networks}} for group TNA models.
#'
#' @export
generate_tna_matrix <- function(nodes_per_group = 5,
                                 group_names = NULL,
                                 n_groups = 3,
                                 edge_prob_range = c(0, 0.3),
                                 self_loops = FALSE,
                                 use_learning_states = TRUE,
                                 learning_categories = NULL,
                                 node_prefix = "N",
                                 seed = NULL,
                                 verbose = TRUE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Determine number of groups
  if (!is.null(group_names)) {
    n_groups <- length(group_names)
  }

  # Handle nodes_per_group: single value, vector matching n_groups, or range
  if (length(nodes_per_group) == 1) {
    nodes_per_group <- rep(nodes_per_group, n_groups)
  } else if (length(nodes_per_group) == 2 && nodes_per_group[1] != nodes_per_group[2]) {
    # Treat as range c(min, max)
    nodes_per_group <- sample(nodes_per_group[1]:nodes_per_group[2], n_groups, replace = TRUE)
  } else if (length(nodes_per_group) > 2 && length(nodes_per_group) != n_groups) {
    # Treat as range, sample for each group
    nodes_per_group <- sample(nodes_per_group, n_groups, replace = TRUE)
  }

  if (length(nodes_per_group) != n_groups) {
    stop("nodes_per_group must be length 1, 2 (range), or match n_groups")
  }

  total_nodes <- sum(nodes_per_group)

  # Generate group names if not provided
  if (is.null(group_names)) {
    group_names <- paste0("Group", 1:n_groups)
  }

  # Determine learning categories per group
  all_cats <- c("metacognitive", "cognitive", "behavioral",
                "social", "motivational", "affective", "group_regulation")

  if (use_learning_states) {
    if (is.null(learning_categories)) {
      # Random category for each group
      learning_categories <- sample(all_cats, n_groups, replace = TRUE)
    } else if (length(learning_categories) == 1) {
      learning_categories <- rep(learning_categories, n_groups)
    } else if (length(learning_categories) != n_groups) {
      stop("learning_categories must be length 1 or match number of groups")
    }
  }

  # Generate node names for each group
  node_types <- list()
  all_nodes <- character(0)

  for (i in seq_len(n_groups)) {
    n_nodes <- nodes_per_group[i]

    if (use_learning_states) {
      node_names <- get_learning_states(learning_categories[i], n = n_nodes)
      # Ensure unique names across groups by adding suffix if needed
      while (any(node_names %in% all_nodes)) {
        duplicates <- node_names %in% all_nodes
        node_names[duplicates] <- paste0(node_names[duplicates], i)
      }
    } else {
      start_idx <- length(all_nodes) + 1
      end_idx <- start_idx + n_nodes - 1
      node_names <- paste0(node_prefix, start_idx:end_idx)
    }

    node_types[[group_names[i]]] <- node_names
    all_nodes <- c(all_nodes, node_names)
  }

  if (verbose) {
    message(sprintf("Generating TNA matrix: %d groups, %d total nodes",
                    n_groups, total_nodes))
    for (i in seq_len(n_groups)) {
      message(sprintf("  %s: %s", group_names[i],
                      paste(node_types[[group_names[i]]], collapse = ", ")))
    }
  }

  # Generate transition matrix
  m <- matrix(
    runif(total_nodes^2, edge_prob_range[1], edge_prob_range[2]),
    total_nodes, total_nodes
  )

  if (!self_loops) {
    diag(m) <- 0
  }

  colnames(m) <- rownames(m) <- all_nodes

  if (verbose) {
    message("TNA matrix generated successfully.")
  }

  return(list(
    matrix = m,
    node_types = node_types
  ))
}
