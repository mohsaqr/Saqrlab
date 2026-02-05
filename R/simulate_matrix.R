#' Simulate Network Matrix
#'
#' @description
#' Generate a simple transition matrix for TNA/Markov network analysis.
#' Each call randomly selects a learning category for node names.
#'
#' @param n_nodes Integer. Number of nodes. Default: 9.
#' @param matrix_type Character. Type of matrix to generate:
#'   \itemize{
#'     \item \code{"transition"}: Directed, row-normalized (rows sum to 1)
#'     \item \code{"frequency"}: Directed, integer counts
#'     \item \code{"co-occurrence"}: Symmetric
#'     \item \code{"adjacency"}: Binary or weighted edges
#'   }
#'   Default: "transition".
#' @param edge_prob Numeric. Probability of edges existing. Default: 0.3.
#' @param weighted Logical. If TRUE, generates weighted edges. Default: TRUE.
#' @param weight_range Numeric vector of length 2. Range for edge weights.
#'   Default: c(0, 1).
#' @param directed Logical. If TRUE, generates directed matrix. Default: TRUE.
#' @param allow_self_loops Logical. If TRUE, allows diagonal entries. Default: FALSE.
#' @param names Character vector or NULL. Custom node names. If NULL, randomly
#'   selects learning states from a random category. Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @return A numeric matrix with row and column names set to learning states.
#'
#' @details
#' This function generates a simple network matrix. Each call randomly picks
#' one learning category (metacognitive, cognitive, behavioral, social,
#' motivational, affective, or group_regulation) and uses verbs from that
#' category as node names.
#'
#' For matrices with multiple node types, use \code{\link{simulate_htna}}.
#'
#' @examples
#' # Simple 9-node transition matrix
#' mat <- simulate_matrix(seed = 42)
#' mat
#' rowSums(mat)  # Rows sum to 1
#'
#' # Frequency matrix
#' mat <- simulate_matrix(n_nodes = 5, matrix_type = "frequency", seed = 42)
#'
#' # Co-occurrence matrix (symmetric)
#' mat <- simulate_matrix(n_nodes = 6, matrix_type = "co-occurrence", seed = 42)
#' isSymmetric(mat)  # TRUE
#'
#' # Custom names
#' mat <- simulate_matrix(n_nodes = 4, names = c("A", "B", "C", "D"), seed = 42)
#'
#' @seealso \code{\link{simulate_htna}} for multi-type matrices
#'
#' @export
simulate_matrix <- function(n_nodes = 9,
                            matrix_type = c("transition", "frequency",
                                            "co-occurrence", "adjacency"),
                            edge_prob = 0.3,
                            weighted = TRUE,
                            weight_range = c(0, 1),
                            directed = TRUE,
                            allow_self_loops = FALSE,
                            names = NULL,
                            seed = NULL) {
  # Match matrix_type argument
  matrix_type <- match.arg(matrix_type)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Auto-adjust directed based on matrix_type
  if (matrix_type == "co-occurrence") {
    directed <- FALSE
  }

  # Auto-adjust weight_range for frequency
  if (matrix_type == "frequency" && identical(weight_range, c(0, 1))) {
    weight_range <- c(1, 100)
  }

  # --- Generate node names ---
  if (!is.null(names)) {
    if (length(names) < n_nodes) {
      stop("names must have at least n_nodes elements.")
    }
    node_names <- names[1:n_nodes]
  } else {
    # Randomly pick one learning category
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    category <- sample(all_cats, 1)
    available_states <- get_learning_states(categories = category)
    if (n_nodes > length(available_states)) {
      extra_needed <- n_nodes - length(available_states)
      extra_names <- paste0("State", seq_len(extra_needed))
      node_names <- c(sample(available_states), extra_names)
    } else {
      node_names <- sample(available_states, n_nodes)
    }
  }

  # --- Initialize matrix ---
  mat <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  rownames(mat) <- node_names
  colnames(mat) <- node_names

  # --- Generate edges ---
  for (i in 1:n_nodes) {
    j_start <- if (directed) 1 else i
    for (j in j_start:n_nodes) {
      if (!allow_self_loops && i == j) next

      if (runif(1) < edge_prob) {
        if (weighted) {
          weight <- runif(1, min = weight_range[1], max = weight_range[2])
          mat[i, j] <- round(weight, 4)
          if (!directed && i != j) {
            mat[j, i] <- mat[i, j]
          }
        } else {
          mat[i, j] <- 1
          if (!directed && i != j) {
            mat[j, i] <- 1
          }
        }
      }
    }
  }

  # For directed, also sample reverse edges
  if (directed) {
    for (i in 1:(n_nodes - 1)) {
      for (j in (i + 1):n_nodes) {
        if (runif(1) < edge_prob) {
          if (weighted) {
            mat[j, i] <- round(runif(1, weight_range[1], weight_range[2]), 4)
          } else {
            mat[j, i] <- 1
          }
        }
      }
    }
  }

  # --- Post-process based on matrix_type ---
  if (matrix_type == "frequency") {
    mat <- round(mat * 100)
    mat <- matrix(as.integer(mat), nrow = n_nodes,
                  dimnames = list(node_names, node_names))
  } else if (matrix_type == "transition") {
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- 1
    mat <- mat / row_sums
    mat <- round(mat, 4)
  } else if (matrix_type == "co-occurrence") {
    mat <- (mat + t(mat)) / 2
    mat <- round(mat, 4)
  }

  return(mat)
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
#' @importFrom stats setNames
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
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Default references
  default_type_names <- c("Metacognitive", "Cognitive", "Behavioral",
                          "Social", "Motivational")
  default_cats <- c("metacognitive", "cognitive", "behavioral",
                    "social", "motivational")

  # Infer n_types from custom type_names
  if (!identical(type_names, default_type_names)) {
    n_types <- length(type_names)
  }

  # Adjust type_names if n_types changed
  if (n_types != 5 && identical(type_names, default_type_names)) {
    all_types <- c("Metacognitive", "Cognitive", "Behavioral",
                   "Social", "Motivational", "Affective", "GroupRegulation")
    type_names <- all_types[1:min(n_types, length(all_types))]
    if (n_types > length(all_types)) {
      type_names <- c(type_names, paste0("Type", (length(all_types) + 1):n_types))
    }
  }

  # Adjust categories to match n_types
  if (n_types != 5 && identical(categories, default_cats)) {
    all_cats <- c("metacognitive", "cognitive", "behavioral",
                  "social", "motivational", "affective", "group_regulation")
    categories <- all_cats[1:min(n_types, length(all_cats))]
    if (n_types > length(all_cats)) {
      categories <- c(categories, rep("all", n_types - length(all_cats)))
    }
  } else if (length(categories) != n_types) {
    categories <- rep_len(categories, n_types)
  }

  # Total nodes
  total_nodes <- n_nodes * n_types
  n_nodes_per_type <- rep(n_nodes, n_types)

  # Generate node names from each category
  node_names <- character(0)
  node_types_list <- list()

  for (i in seq_len(n_types)) {
    states <- get_learning_states(categories = categories[i])
    if (n_nodes > length(states)) {
      extra <- paste0(type_names[i], "_", seq_len(n_nodes - length(states)))
      type_nodes <- c(sample(states), extra)
    } else {
      type_nodes <- sample(states, n_nodes)
    }
    node_names <- c(node_names, type_nodes)
    node_types_list[[type_names[i]]] <- type_nodes
  }

  # Build type indices for edge probability
  node_type_idx <- rep(seq_len(n_types), each = n_nodes)

  # Initialize matrix
  mat <- matrix(0, nrow = total_nodes, ncol = total_nodes)
  rownames(mat) <- node_names
  colnames(mat) <- node_names

  # Generate edges with within/between probabilities
  for (i in 1:total_nodes) {
    for (j in 1:total_nodes) {
      if (!allow_self_loops && i == j) next

      # Determine edge probability
      if (node_type_idx[i] == node_type_idx[j]) {
        prob <- within_prob
      } else {
        prob <- between_prob
      }

      if (runif(1) < prob) {
        mat[i, j] <- round(runif(1, weight_range[1], weight_range[2]), 4)
      }
    }
  }

  # Normalize rows to sum to 1 (transition matrix)
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  mat <- mat / row_sums
  mat <- round(mat, 4)

  return(list(
    matrix = mat,
    node_types = node_types_list,
    type_names = type_names,
    n_nodes_per_type = setNames(n_nodes_per_type, type_names)
  ))
}


#' @rdname simulate_htna
#' @export
simulate_mlna <- simulate_htna


#' @rdname simulate_htna
#' @export
simulate_mtna <- simulate_htna
