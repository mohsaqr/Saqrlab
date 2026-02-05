#' Simulate Network Matrix
#'
#' @description
#' Generate a simulated matrix for network analysis. Supports different matrix
#' types (transition, frequency, co-occurrence, adjacency) and multiple node
#' types (communities/groups) with configurable connection probabilities.
#'
#' @param n_nodes Integer or integer vector. If single value, total number of nodes.
#'   If vector, number of nodes per type. Default: 9.
#' @param n_types Integer. Number of node types/groups. Ignored if `n_nodes` is a
#'   vector. Default: 1.
#' @param matrix_type Character. Type of matrix to generate:
#'   \itemize{
#'     \item \code{"transition"}: Directed, row-normalized (rows sum to 1), like Markov/TNA
#'     \item \code{"frequency"}: Directed, integer counts of transitions
#'     \item \code{"co-occurrence"}: Symmetric, counts of co-occurrences
#'     \item \code{"adjacency"}: Binary or weighted edges (default behavior)
#'   }
#'   Default: "transition".
#' @param within_prob Numeric or numeric vector. Probability of edges within each
#'   type. If single value, applies to all types. Default: 0.3.
#' @param between_prob Numeric or matrix. Probability of edges between different
#'   types. If single value, applies to all between-type pairs. If matrix, should
#'   be symmetric with dimensions n_types x n_types. Default: 0.1.
#' @param weighted Logical. If TRUE, generates weighted edges. If FALSE, binary
#'   (0/1) matrix. Ignored for "transition" and "frequency" types. Default: TRUE.
#' @param weight_range Numeric vector of length 2. Range for edge weights.
#'   For "transition": c(0, 1). For "frequency": c(1, 100). For others: c(0, 1).
#'   Default: NULL (auto-set based on matrix_type).
#' @param directed Logical. If TRUE, generates directed (asymmetric) matrix.
#'   Auto-set to TRUE for "transition"/"frequency", FALSE for "co-occurrence".
#'   Default: NULL (auto-set based on matrix_type).
#' @param allow_self_loops Logical. If TRUE, allows diagonal entries (self-loops).
#'   Default: FALSE.
#' @param normalize Logical. If TRUE for "transition" type, normalizes rows to
#'   sum to 1. Default: TRUE.
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
#'   \item{matrix_type}{The type of matrix generated.}
#' }
#'
#' @details
#' The function implements a stochastic block model approach:
#' \itemize{
#'   \item Nodes are assigned to types/groups
#'   \item Edges within the same type occur with probability `within_prob`
#'   \item Edges between different types occur with probability `between_prob`
#'   \item Edge weights depend on `matrix_type`
#' }
#'
#' **Matrix Types**:
#' \itemize{
#'   \item \strong{transition}: Rows sum to 1 (Markov chain style). Suitable for
#'     TNA, sequence analysis, and state transition modeling.
#'   \item \strong{frequency}: Integer counts representing transition frequencies.
#'     Useful for raw count data before normalization.
#'   \item \strong{co-occurrence}: Symmetric matrix of co-occurrence counts.
#'     Useful for undirected association networks.
#'   \item \strong{adjacency}: Standard adjacency matrix (binary or weighted).
#' }
#'
#' @examples
#' # Default: 9-node transition matrix (TNA style)
#' result <- simulate_matrix(seed = 42)
#' result$matrix
#' rowSums(result$matrix)  # Rows sum to ~1
#'
#' # Frequency matrix (integer counts)
#' result <- simulate_matrix(
#'   n_nodes = 5,
#'   matrix_type = "frequency",
#'   seed = 42
#' )
#' result$matrix  # Integer values
#'
#' # Co-occurrence matrix (symmetric)
#' result <- simulate_matrix(
#'   n_nodes = 6,
#'   matrix_type = "co-occurrence",
#'   seed = 42
#' )
#' isSymmetric(result$matrix)  # TRUE
#'
#' # Adjacency matrix (binary)
#' result <- simulate_matrix(
#'   n_nodes = 8,
#'   matrix_type = "adjacency",
#'   weighted = FALSE,
#'   seed = 42
#' )
#'
#' # Multi-type network
#' result <- simulate_matrix(
#'   n_nodes = c(3, 3, 3),
#'   type_names = c("Meta", "Cog", "Behav"),
#'   within_prob = 0.5,
#'   between_prob = 0.2,
#'   seed = 42
#' )
#'
#' @seealso \code{\link{simulate_htna}}, \code{\link{simulate_edge_list}}
#'
#' @importFrom stats setNames
#' @export
simulate_matrix <- function(n_nodes = 9,
                            n_types = 1,
                            matrix_type = c("transition", "frequency",
                                            "co-occurrence", "adjacency"),
                            within_prob = 0.3,
                            between_prob = 0.1,
                            weighted = TRUE,
                            weight_range = NULL,
                            directed = NULL,
                            allow_self_loops = FALSE,
                            normalize = TRUE,
                            use_learning_states = TRUE,
                            categories = "all",
                            names = NULL,
                            type_names = NULL,
                            seed = NULL) {
  # Match matrix_type argument
 matrix_type <- match.arg(matrix_type)

  # Set defaults based on matrix_type
  if (is.null(directed)) {
    directed <- matrix_type %in% c("transition", "frequency", "adjacency")
  }
  if (is.null(weight_range)) {
    weight_range <- switch(matrix_type,
      "transition" = c(0, 1),
      "frequency" = c(1, 100),
      "co-occurrence" = c(1, 50),
      "adjacency" = c(0, 1)
    )
  }
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

  # --- Post-process based on matrix_type ---
  if (matrix_type == "frequency") {
    # Convert to integers
    adj_matrix <- round(adj_matrix * 100)
    adj_matrix <- matrix(as.integer(adj_matrix), nrow = total_nodes,
                         dimnames = list(node_names, node_names))
  } else if (matrix_type == "transition" && normalize) {
    # Normalize rows to sum to 1
    row_sums <- rowSums(adj_matrix)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    adj_matrix <- adj_matrix / row_sums
    adj_matrix <- round(adj_matrix, 4)
  } else if (matrix_type == "co-occurrence") {
    # Ensure symmetric
    adj_matrix <- (adj_matrix + t(adj_matrix)) / 2
    adj_matrix <- round(adj_matrix, 4)
  }

  # --- Return results ---
  return(list(
    matrix = adj_matrix,
    node_types = node_types,
    type_names = type_names,
    n_nodes_per_type = setNames(n_nodes_per_type, type_names),
    matrix_type = matrix_type
  ))
}


#' Simulate HTNA/MLNA/MTNA Matrix with Node Types
#'
#' @description
#' Generate a transition matrix with multiple node types for hierarchical (HTNA),
#' multilevel (MLNA), or multi-type (MTNA) network analysis. By default creates
#' a 25-node matrix (5 nodes x 5 types) using learning category names.
#'
#' `simulate_htna()`, `simulate_mlna()`, and `simulate_mtna()` are aliases that
#' produce identical output - use whichever name fits your analysis context.
#'
#' @param n_nodes Integer. Number of nodes per type. Default: 5.
#' @param n_types Integer. Number of node types. Default: 5.
#' @param type_names Character vector. Names for node types. Default uses
#'   learning categories: "Metacognitive", "Cognitive", "Behavioral",
#'   "Social", "Motivational".
#' @param within_prob Numeric. Probability of edges within each type. Default: 0.4.
#' @param between_prob Numeric. Probability of edges between types. Default: 0.15.
#' @param weight_range Numeric vector of length 2. Range for edge weights.
#'   Default: c(0, 1).
#' @param allow_self_loops Logical. Allow diagonal entries. Default: FALSE.
#' @param categories Character vector. Learning state categories for node names.
#'   Can be single (same for all types) or one per type. Default: matches type_names.
#' @param seed Integer or NULL. Random seed. Default: NULL.
#'
#' @return A list containing:
#' \describe{
#'   \item{matrix}{Numeric transition matrix (rows sum to 1).}
#'   \item{node_types}{Named list mapping type names to node names (for plot_htna/plot_mlna).}
#'   \item{type_names}{Character vector of type names.}
#'   \item{n_nodes_per_type}{Integer vector of node counts per type.}
#' }
#'
#' @examples
#' # Default: 5 types x 5 nodes = 25 node matrix
#' net <- simulate_htna(seed = 42)
#' net$matrix
#' net$node_types  # List format for plot_htna/plot_mlna
#'
#' # Custom number of types
#' net <- simulate_htna(n_nodes = 4, n_types = 3, seed = 42)
#'
#' # Custom type names
#' net <- simulate_mlna(
#'   n_nodes = 6,
#'   type_names = c("Macro", "Meso", "Micro"),
#'   seed = 42
#' )
#'
#' # Use with tna package:
#' # plot_htna(net$matrix, net$node_types, layout = "polygon")
#' # plot_mlna(net$matrix, layers = net$node_types)
#'
#' @seealso \code{\link{simulate_matrix}}, \code{\link{generate_tna_matrix}}
#'
#' @name simulate_htna
#' @rdname simulate_htna
#' @export
simulate_htna <- function(n_nodes = 5,
                          n_types = 5,
                          type_names = c("Metacognitive", "Cognitive",
                                         "Behavioral", "Social", "Motivational"),
                          within_prob = 0.4,
                          between_prob = 0.15,
                          weight_range = c(0, 1),
                          allow_self_loops = FALSE,
                          categories = c("metacognitive", "cognitive",
                                         "behavioral", "social", "motivational"),
                          seed = NULL) {
  # Default type/category names for reference
  default_type_names <- c("Metacognitive", "Cognitive", "Behavioral",
                          "Social", "Motivational")
  default_cats <- c("metacognitive", "cognitive", "behavioral",
                    "social", "motivational")

  # If custom type_names provided, infer n_types from it
  if (!identical(type_names, default_type_names)) {
    n_types <- length(type_names)
  }

  # Adjust type_names if n_types differs from default but type_names wasn't customized
  if (n_types != 5 && identical(type_names, default_type_names)) {
    all_types <- c("Metacognitive", "Cognitive", "Behavioral",
                   "Social", "Motivational", "Affective", "GroupRegulation")
    type_names <- all_types[1:min(n_types, length(all_types))]
    if (n_types > length(all_types)) {
      type_names <- c(type_names, paste0("Type", (length(all_types) + 1):n_types))
    }
  }

  # Adjust categories if n_types differs and categories wasn't customized
  if (n_types != 5 && identical(categories, default_cats)) {
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    categories <- all_cats[1:min(n_types, length(all_cats))]
    if (n_types > length(all_cats)) {
      categories <- c(categories, rep("all", n_types - length(all_cats)))
    }
  } else if (length(categories) != n_types) {
    # If categories provided but wrong length, recycle or truncate
    if (length(categories) == 1) {
      categories <- rep(categories, n_types)
    } else {
      categories <- rep_len(categories, n_types)
    }
  }

  # Build n_nodes vector (one per type)
  n_nodes_vec <- rep(n_nodes, n_types)

  # Call simulate_matrix
  result <- simulate_matrix(
    n_nodes = n_nodes_vec,
    n_types = n_types,
    matrix_type = "transition",
    within_prob = within_prob,
    between_prob = between_prob,
    weighted = TRUE,
    weight_range = weight_range,
    directed = TRUE,
    allow_self_loops = allow_self_loops,
    normalize = TRUE,
    use_learning_states = TRUE,
    categories = categories,
    type_names = type_names,
    seed = seed
  )

  # Convert node_types from named vector to list format (for HTNA/MLNA compatibility)
  node_types_list <- list()
  for (type in result$type_names) {
    node_types_list[[type]] <- names(result$node_types[result$node_types == type])
  }

  return(list(
    matrix = result$matrix,
    node_types = node_types_list,
    type_names = result$type_names,
    n_nodes_per_type = result$n_nodes_per_type
  ))
}


#' @rdname simulate_htna
#' @export
simulate_mlna <- simulate_htna


#' @rdname simulate_htna
#' @export
simulate_mtna <- simulate_htna
