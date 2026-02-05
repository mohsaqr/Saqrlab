#' Simulate Network Matrix with Node Types
#'
#' @description
#' Generate a simulated adjacency or weight matrix for network analysis.
#' Supports single or multiple node types (communities/groups) with configurable
#' within-group and between-group connection probabilities.
#'
#' @param n_nodes Integer or integer vector. If single value, total number of nodes.
#'   If vector, number of nodes per type. Default: 20.
#' @param n_types Integer. Number of node types/groups. Ignored if `n_nodes` is a
#'   vector. Default: 1.
#' @param within_prob Numeric or numeric vector. Probability of edges within each
#'   type. If single value, applies to all types. Default: 0.3.
#' @param between_prob Numeric or matrix. Probability of edges between different
#'   types. If single value, applies to all between-type pairs. If matrix, should
#'   be symmetric with dimensions n_types x n_types. Default: 0.1.
#' @param weighted Logical. If TRUE, generates weighted edges. If FALSE, binary
#'   (0/1) matrix. Default: TRUE.
#' @param weight_range Numeric vector of length 2. Range for edge weights when
#'   `weighted = TRUE`. Default: c(0.1, 1.0).
#' @param directed Logical. If TRUE, generates directed (asymmetric) matrix.
#'   Default: FALSE.
#' @param allow_self_loops Logical. If TRUE, allows diagonal entries (self-loops).
#'   Default: FALSE.
#' @param use_learning_states Logical. If TRUE, uses learning action verbs as
#'   node names. Default: TRUE.
#' @param categories Character vector. Learning state categories to use when
#'   `use_learning_states = TRUE`. Options: "metacognitive", "cognitive",
#'   "behavioral", "social", "motivational", "affective", "group_regulation",
#'   or "all". Default: "all".
#' @param names Character vector or NULL. Custom node names. Overrides
#'   `use_learning_states` if provided. If NULL and `use_learning_states = FALSE`,
#'   uses names from `GLOBAL_NAMES`. Default: NULL.
#' @param type_names Character vector or NULL. Names for node types. If NULL,
#'   uses "Type1", "Type2", etc. Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @return A list containing:
#' \describe{
#'   \item{matrix}{Numeric matrix with row and column names set to node names.}
#'   \item{node_types}{Named character vector mapping node names to their types.}
#'   \item{type_names}{Character vector of type names.}
#'   \item{n_nodes_per_type}{Integer vector of node counts per type.}
#' }
#'
#' @details
#' The function implements a stochastic block model approach:
#' \itemize{
#'   \item Nodes are assigned to types/groups
#'   \item Edges within the same type occur with probability `within_prob`
#'   \item Edges between different types occur with probability `between_prob`
#'   \item Edge weights (if enabled) are uniformly distributed in `weight_range`
#' }
#'
#' For undirected networks, the matrix is symmetric. For directed networks,
#' each potential edge is sampled independently.
#'
#' @examples
#' # Simple matrix with learning states (default)
#' result <- simulate_matrix(n_nodes = 10, seed = 42)
#' result$matrix
#'
#' # Matrix with specific learning categories
#' result <- simulate_matrix(
#'   n_nodes = 8,
#'   categories = c("metacognitive", "cognitive"),
#'   seed = 42
#' )
#' colnames(result$matrix)
#'
#' # Two-type network with learning states
#' result <- simulate_matrix(
#'   n_nodes = c(6, 4),
#'   type_names = c("Regulation", "Cognition"),
#'   within_prob = c(0.5, 0.6),
#'   between_prob = 0.2,
#'   categories = "group_regulation",
#'   seed = 42
#' )
#' result$matrix
#' table(result$node_types)
#'
#' # Three communities with different densities
#' result <- simulate_matrix(
#'   n_nodes = c(5, 5, 5),
#'   type_names = c("Group_A", "Group_B", "Group_C"),
#'   within_prob = c(0.5, 0.4, 0.3),
#'   between_prob = 0.05,
#'   weighted = FALSE,
#'   seed = 123
#' )
#'
#' # Using global names instead of learning states
#' result <- simulate_matrix(
#'   n_nodes = 10,
#'   use_learning_states = FALSE,
#'   seed = 42
#' )
#' colnames(result$matrix)  # Will show names like "Emma", "Liam", etc.
#'
#' # Directed network
#' result <- simulate_matrix(
#'   n_nodes = 12,
#'   n_types = 2,
#'   directed = TRUE,
#'   seed = 42
#' )
#'
#' @seealso \code{\link{simulate_edge_list}}, \code{\link{GLOBAL_NAMES}}
#'
#' @importFrom stats setNames
#' @export
simulate_matrix <- function(n_nodes = 20,
                            n_types = 1,
                            within_prob = 0.3,
                            between_prob = 0.1,
                            weighted = TRUE,
                            weight_range = c(0.1, 1.0),
                            directed = FALSE,
                            allow_self_loops = FALSE,
                            use_learning_states = TRUE,
                            categories = "all",
                            names = NULL,
                            type_names = NULL,
                            seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # --- Determine node counts per type ---

if (length(n_nodes) > 1) {
    # n_nodes is a vector specifying nodes per type
    n_nodes_per_type <- as.integer(n_nodes)
    n_types <- length(n_nodes_per_type)
    total_nodes <- sum(n_nodes_per_type)
  } else {
    # Single n_nodes value, distribute across types
    total_nodes <- as.integer(n_nodes)
    if (n_types == 1) {
      n_nodes_per_type <- total_nodes
    } else {
      # Distribute nodes as evenly as possible
      base_count <- total_nodes %/% n_types
      remainder <- total_nodes %% n_types
      n_nodes_per_type <- rep(base_count, n_types)
      if (remainder > 0) {
        n_nodes_per_type[1:remainder] <- n_nodes_per_type[1:remainder] + 1
      }
    }
  }

  # --- Input Validation ---
  stopifnot(
    all(n_nodes_per_type >= 1),
    is.numeric(within_prob), all(within_prob >= 0), all(within_prob <= 1),
    is.numeric(between_prob), all(between_prob >= 0), all(between_prob <= 1),
    is.logical(weighted), length(weighted) == 1,
    is.numeric(weight_range), length(weight_range) == 2,
    weight_range[1] < weight_range[2],
    is.logical(directed), length(directed) == 1,
    is.logical(allow_self_loops), length(allow_self_loops) == 1
  )

  # --- Expand within_prob if needed ---
  if (length(within_prob) == 1) {
    within_prob <- rep(within_prob, n_types)
  } else if (length(within_prob) != n_types) {
    stop("within_prob must be length 1 or equal to number of types.")
  }

  # --- Build probability matrix ---
  if (is.matrix(between_prob)) {
    if (nrow(between_prob) != n_types || ncol(between_prob) != n_types) {
      stop("between_prob matrix must have dimensions n_types x n_types.")
    }
    prob_matrix <- between_prob
    # Override diagonal with within_prob
    diag(prob_matrix) <- within_prob
  } else {
    # Single between_prob value
    prob_matrix <- matrix(between_prob, nrow = n_types, ncol = n_types)
    diag(prob_matrix) <- within_prob
  }

  # --- Generate type names ---
  if (is.null(type_names)) {
    type_names <- paste0("Type", 1:n_types)
  } else if (length(type_names) != n_types) {
    stop("type_names must have length equal to number of types.")
  }

  # --- Generate node names ---
  if (!is.null(names)) {
    # Custom names provided
    if (length(names) < total_nodes) {
      stop("Provided names vector must have at least total_nodes elements.")
    }
    node_names <- names[1:total_nodes]
  } else if (use_learning_states) {
    # Use learning states
    available_states <- get_learning_states(categories = categories)
    if (total_nodes > length(available_states)) {
      # Need more names than available states
      extra_needed <- total_nodes - length(available_states)
      extra_names <- paste0("State", seq_len(extra_needed))
      node_names <- c(sample(available_states), extra_names)
    } else {
      node_names <- sample(available_states, total_nodes)
    }
  } else {
    # Use global names
    if (total_nodes > length(GLOBAL_NAMES)) {
      base_names <- GLOBAL_NAMES
      extra_needed <- total_nodes - length(GLOBAL_NAMES)
      extra_names <- paste0("Node", seq_len(extra_needed))
      node_names <- c(sample(base_names), extra_names)
    } else {
      node_names <- sample(GLOBAL_NAMES, total_nodes)
    }
  }

  # --- Assign nodes to types ---
  node_type_indices <- rep(1:n_types, times = n_nodes_per_type)
  node_types <- setNames(type_names[node_type_indices], node_names)

  # --- Initialize matrix ---
  adj_matrix <- matrix(0, nrow = total_nodes, ncol = total_nodes)
  rownames(adj_matrix) <- node_names
  colnames(adj_matrix) <- node_names

  # --- Generate edges ---
  for (i in 1:total_nodes) {
    j_start <- if (directed) 1 else i
    for (j in j_start:total_nodes) {
      # Skip self-loops if not allowed
      if (!allow_self_loops && i == j) next

      # Get types
      type_i <- node_type_indices[i]
      type_j <- node_type_indices[j]

      # Get probability for this pair
      edge_prob <- prob_matrix[type_i, type_j]

      # Sample edge
      if (runif(1) < edge_prob) {
        if (weighted) {
          weight <- runif(1, min = weight_range[1], max = weight_range[2])
          adj_matrix[i, j] <- round(weight, 4)
          if (!directed && i != j) {
            adj_matrix[j, i] <- adj_matrix[i, j]
          }
        } else {
          adj_matrix[i, j] <- 1
          if (!directed && i != j) {
            adj_matrix[j, i] <- 1
          }
        }
      }
    }
  }

  # --- For directed networks, also sample j->i edges ---
  if (directed) {
    for (i in 1:total_nodes) {
      for (j in 1:total_nodes) {
        if (i == j && !allow_self_loops) next
        if (i >= j) next  # Already processed in the loop above for i->j

        type_i <- node_type_indices[i]
        type_j <- node_type_indices[j]
        edge_prob <- prob_matrix[type_j, type_i]  # Note: reversed for j->i

        if (runif(1) < edge_prob) {
          if (weighted) {
            weight <- runif(1, min = weight_range[1], max = weight_range[2])
            adj_matrix[j, i] <- round(weight, 4)
          } else {
            adj_matrix[j, i] <- 1
          }
        }
      }
    }
  }

  # --- Return results ---
  return(list(
    matrix = adj_matrix,
    node_types = node_types,
    type_names = type_names,
    n_nodes_per_type = setNames(n_nodes_per_type, type_names)
  ))
}
