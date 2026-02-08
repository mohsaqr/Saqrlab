#' Simulate Social Network Edge List
#'
#' @description
#' Generate a simulated edge list for social network analysis. Creates random
#' connections between nodes with weights and class assignments.
#'
#' @param n_nodes Integer. Number of nodes (people) in the network. Default: 20.
#' @param n_edges Integer or NULL. Number of edges to generate. If NULL,
#'   calculated as `n_nodes * edge_density`. Default: NULL.
#' @param edge_density Numeric. Average number of edges per node when `n_edges`
#'   is NULL. Default: 3.
#' @param n_classes Integer. Number of classes/groups (2-10). Default: 3.
#' @param directed Logical. Whether edges are directed. Default: TRUE.
#' @param allow_self_loops Logical. Whether to allow self-connections. Default: FALSE.
#' @param weight_range Numeric vector of length 2. Range for edge weights.
#'   Default: c(0.1, 1.0).
#' @param names Character vector or NULL. Custom node names. If NULL, uses
#'   names from `GLOBAL_NAMES`. Default: NULL.
#' @param class_probs Numeric vector or NULL. Probability of each class. Must
#'   sum to 1 and have length `n_classes`. If NULL, uniform distribution.
#'   Default: NULL.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{source}{Character. Name of the source node.}
#'   \item{target}{Character. Name of the target node.}
#'   \item{weight}{Numeric. Edge weight in the specified range.}
#'   \item{class}{Integer. Class assignment (1 to n_classes).}
#' }
#'
#' @details
#' The function generates a random social network edge list with:
#' \itemize{
#'   \item Nodes named using diverse global names (or custom names)
#'   \item Random edges between nodes
#'   \item Weights uniformly distributed in the specified range
#'   \item Class assignments based on specified probabilities
#' }
#'
#' For undirected networks, each edge appears once (no duplicate A-B, B-A pairs).
#'
#' @examples
#' # Basic usage with defaults
#' edges <- simulate_edge_list(seed = 42)
#' head(edges)
#'
#' # Larger network with 5 classes
#' edges <- simulate_edge_list(
#'   n_nodes = 50,
#'   n_classes = 5,
#'   edge_density = 4,
#'   seed = 123
#' )
#'
#' # Custom names and specific number of edges
#' edges <- simulate_edge_list(
#'   n_nodes = 10,
#'   n_edges = 30,
#'   names = c("Alice", "Bob", "Carol", "Dave", "Eve",
#'             "Frank", "Grace", "Hank", "Ivy", "Jack"),
#'   n_classes = 2,
#'   seed = 42
#' )
#'
#' # Undirected network with unequal class distribution
#' edges <- simulate_edge_list(
#'   n_nodes = 30,
#'   directed = FALSE,
#'   n_classes = 3,
#'   class_probs = c(0.5, 0.3, 0.2),
#'   seed = 42
#' )
#'
#' # View class distribution
#' table(edges$class)
#'
#' @seealso \code{\link{GLOBAL_NAMES}}, \code{\link{get_global_names}}
#'
#' @export
simulate_edge_list <- function(n_nodes = 20,
                                n_edges = NULL,
                                edge_density = 3,
                                n_classes = 3,
                                directed = TRUE,
                                allow_self_loops = FALSE,
                                weight_range = c(0.1, 1.0),
                                names = NULL,
                                class_probs = NULL,
                                seed = NULL) {
# Set seed if provided
if (!is.null(seed)) set.seed(seed)

# --- Input Validation ---
stopifnot(
    is.numeric(n_nodes), length(n_nodes) == 1, n_nodes >= 2,
    is.numeric(edge_density), length(edge_density) == 1, edge_density > 0,
    is.numeric(n_classes), length(n_classes) == 1, n_classes >= 2, n_classes <= 10,
    is.logical(directed), length(directed) == 1,
    is.logical(allow_self_loops), length(allow_self_loops) == 1,
    is.numeric(weight_range), length(weight_range) == 2,
    weight_range[1] < weight_range[2]
)

# Validate class_probs
if (!is.null(class_probs)) {
    if (length(class_probs) != n_classes) {
    stop("class_probs must have length equal to n_classes.")
    }
    if (abs(sum(class_probs) - 1) > 1e-6) {
    stop("class_probs must sum to 1.")
    }
}

# --- Determine node names ---
if (is.null(names)) {
    node_names <- get_global_names(n_nodes)
} else {
    if (length(names) < n_nodes) {
    stop("Provided names vector must have at least n_nodes elements.")
    }
    node_names <- names[1:n_nodes]
}

# --- Calculate number of edges ---
if (is.null(n_edges)) {
    n_edges <- round(n_nodes * edge_density)
}

# Calculate maximum possible edges
if (directed) {
    max_edges <- if (allow_self_loops) n_nodes^2 else n_nodes * (n_nodes - 1)
} else {
    max_edges <- if (allow_self_loops) {
    n_nodes * (n_nodes + 1) / 2
    } else {
    n_nodes * (n_nodes - 1) / 2
    }
}

if (n_edges > max_edges) {
    warning(sprintf(
    "Requested %d edges exceeds maximum possible (%d). Using maximum.",
    n_edges, max_edges
    ))
    n_edges <- max_edges
}

# --- Generate edges ---
edges_set <- list()
edge_count <- 0
attempts <- 0
max_attempts <- n_edges * 100 # Prevent infinite loops

while (edge_count < n_edges && attempts < max_attempts) {
    attempts <- attempts + 1

    # Random source and target
    source_idx <- sample.int(n_nodes, 1)
    target_idx <- sample.int(n_nodes, 1)

    # Skip self-loops if not allowed
    if (!allow_self_loops && source_idx == target_idx) next

    # For undirected, normalize edge order
    if (!directed && source_idx > target_idx) {
    temp <- source_idx
    source_idx <- target_idx
    target_idx <- temp
    }

    # Create edge key
    edge_key <- paste(source_idx, target_idx, sep = "-")

    # Skip if edge already exists
    if (edge_key %in% names(edges_set)) next

    # Add edge
    edges_set[[edge_key]] <- c(source_idx, target_idx)
    edge_count <- edge_count + 1
}

if (edge_count < n_edges) {
    warning(sprintf(
    "Could only generate %d of %d requested edges.",
    edge_count, n_edges
    ))
}

# --- Build edge list data frame ---
if (edge_count == 0) {
    return(data.frame(
    source = character(),
    target = character(),
    weight = numeric(),
    class = integer(),
    stringsAsFactors = FALSE
    ))
}

edge_matrix <- do.call(rbind, edges_set)

# Generate weights
weights <- runif(edge_count, min = weight_range[1], max = weight_range[2])

# Generate classes
if (is.null(class_probs)) {
    classes <- sample.int(n_classes, edge_count, replace = TRUE)
} else {
    classes <- sample.int(n_classes, edge_count, replace = TRUE, prob = class_probs)
}

# Create data frame
edge_list <- data.frame(
    source = node_names[edge_matrix[, 1]],
    target = node_names[edge_matrix[, 2]],
    weight = round(weights, 4),
    class = classes,
    stringsAsFactors = FALSE
)

# Sort by source, then target
edge_list <- edge_list[order(edge_list$source, edge_list$target), ]
rownames(edge_list) <- NULL

return(edge_list)
}
