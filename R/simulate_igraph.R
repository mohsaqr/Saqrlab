#' Simulate igraph Network Object
#'
#' @description
#' Generate an igraph network object using common graph algorithms with
#' realistic node names from human names or learning states.
#'
#' @param n Integer or NULL. Number of nodes. If NULL (default), randomly
#'   selects between 20-50 nodes.
#' @param model Character. Graph generation algorithm:
#'   \itemize{
#'     \item \code{"er"}: Erdos-Renyi random graph
#'     \item \code{"ba"}: Barabasi-Albert scale-free network
#'     \item \code{"ws"}: Watts-Strogatz small-world network
#'     \item \code{"sbm"}: Stochastic Block Model (community structure)
#'     \item \code{"reg"}: Regular graph (fixed degree)
#'     \item \code{"grg"}: Geometric Random Graph (spatial)
#'     \item \code{"ff"}: Forest Fire (growing network)
#'   }
#'   Default: "er".
#' @param name_source Character. Source for node names:
#'   \itemize{
#'     \item \code{"human"}: Culturally diverse human names from GLOBAL_NAMES
#'     \item \code{"states"}: Learning action verbs from LEARNING_STATES
#'   }
#'   Default: "human".
#' @param regions Character vector. Regions to sample human names from (only used
#'   when name_source = "human"). Can be specific regions (e.g., "arab", "east_asia"),
#'   shortcuts (e.g., "europe", "africa", "asia"), or "all". See \code{\link{list_name_regions}}.
#'   Default: "all".
#' @param categories Character vector. Learning state categories (only used when
#'   name_source = "states"). Options: "metacognitive", "cognitive", "behavioral",
#'   "social", "motivational", "affective", "group_regulation", "lms", or "all".
#'   Default: "all".
#' @param names Character vector or NULL. Custom node names. Overrides
#'   name_source if provided. Default: NULL.
#' @param directed Logical. If TRUE, generate directed network. Default: FALSE.
#' @param weighted Logical. If TRUE, add random edge weights. Default: FALSE.
#' @param weights Numeric vector of length 2. Weight range \[min, max\].
#'   Default: c(0.1, 1.0).
#' @param p Numeric. Edge probability for Erdos-Renyi model. Default: 0.1.
#' @param m Integer or NULL. Fixed number of edges for Erdos-Renyi.
#'   Overrides p if provided. Default: NULL.
#' @param power Numeric. Attachment power for Barabasi-Albert. Default: 1.
#' @param m_ba Integer. Edges per new vertex for Barabasi-Albert. Default: 2.
#' @param nei Integer. Neighborhood size for Watts-Strogatz. Default: 2.
#' @param p_rewire Numeric. Rewiring probability for Watts-Strogatz. Default: 0.05.
#' @param blocks Integer. Number of blocks for SBM. Default: 3.
#' @param p_within Numeric. Within-block edge probability for SBM. Default: 0.3.
#' @param p_between Numeric. Between-block edge probability for SBM. Default: 0.05.
#' @param k Integer. Degree for regular graphs. Default: 4.
#' @param radius Numeric. Connection radius for geometric random graph. Default: 0.25.
#' @param fw Numeric. Forward burning probability for Forest Fire. Default: 0.35.
#' @param bw Numeric. Backward burning factor for Forest Fire. Default: 0.32.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: NULL.
#'
#' @return An igraph object with vertex names and optional edge weights.
#'
#' @details
#' This function wraps igraph's graph generation algorithms and adds:
#' \itemize{
#'   \item Meaningful node names (human names or learning states)
#'   \item Optional edge weights
#'   \item Block/community attributes for SBM
#' }
#'
#' The generated networks are suitable for social network analysis,
#' teaching network concepts, or testing TNA methods.
#'
#' @examples
#' library(igraph)
#'
#' # Default: Erdos-Renyi with human names from all regions
#' g <- simulate_igraph(n = 15, seed = 42)
#' V(g)$name  # Diverse names like "Yuki", "Omar", "Priya"
#' plot(g)
#'
#' # Names from specific regions
#' g_arab <- simulate_igraph(n = 15, regions = "arab", seed = 42)
#' V(g_arab)$name  # Arab names
#'
#' g_africa <- simulate_igraph(n = 20, regions = "africa", seed = 42)
#' V(g_africa)$name  # African names from all sub-regions
#'
#' # Use learning state names instead
#' g_states <- simulate_igraph(n = 10, name_source = "states", seed = 42)
#' V(g_states)$name  # Action verbs like "Plan", "Monitor", "Evaluate"
#'
#' # Scale-free network (Barabasi-Albert)
#' g_sf <- simulate_igraph(n = 50, model = "ba", m_ba = 2, seed = 42)
#' degree_distribution(g_sf)
#'
#' # Small-world network (Watts-Strogatz)
#' g_sw <- simulate_igraph(n = 30, model = "ws", nei = 4, p_rewire = 0.1, seed = 42)
#'
#' # Community structure (Stochastic Block Model)
#' g_sbm <- simulate_igraph(n = 30, model = "sbm", blocks = 3,
#'                          p_within = 0.5, p_between = 0.05, seed = 42)
#' V(g_sbm)$block  # Community assignments
#'
#' # Weighted network with custom names
#' g_custom <- simulate_igraph(
#'   n = 10,
#'   model = "er",
#'   weighted = TRUE,
#'   names = c("Alice", "Bob", "Carol", "Dave", "Eve",
#'             "Frank", "Grace", "Heidi", "Ivan", "Judy"),
#'   seed = 42
#' )
#' E(g_custom)$weight
#'
#' @seealso
#' \code{\link{simulate_network}} for statnet network objects,
#' \code{\link{simulate_matrix}} for adjacency matrices,
#' \code{\link{simulate_tna_network}} for fitted TNA models.
#'
#' @importFrom igraph sample_gnp sample_gnm sample_pa sample_smallworld
#' @importFrom igraph sample_sbm sample_k_regular sample_grg sample_forestfire
#' @importFrom igraph V V<- E E<- ecount is_directed
#' @export
simulate_igraph <- function(n = NULL,
                            model = c("er", "ba", "ws", "sbm", "reg", "grg", "ff"),
                            name_source = c("human", "states"),
                            regions = "all",
                            categories = "all",
                            names = NULL,
                            directed = FALSE,
                            weighted = FALSE,
                            weights = c(0.1, 1.0),
                            p = 0.1,
                            m = NULL,
                            power = 1,
                            m_ba = 2,
                            nei = 2,
                            p_rewire = 0.05,
                            blocks = 3,
                            p_within = 0.3,
                            p_between = 0.05,
                            k = 4,
                            radius = 0.25,
                            fw = 0.35,
                            bw = 0.32,
                            seed = NULL) {
  # Match arguments
  model <- match.arg(model)
  name_source <- match.arg(name_source)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Random n if not specified
  if (is.null(n)) {
    n <- sample(20:50, 1)
  }

  # Generate base graph using internal helper
  g <- .generate_graph(
    n = n, model = model, directed = directed, weighted = weighted,
    weights = weights, p = p, m = m, power = power, m_ba = m_ba,
    nei = nei, p_rewire = p_rewire, blocks = blocks, p_within = p_within,
    p_between = p_between, k = k, radius = radius, fw = fw, bw = bw,
    name_source = name_source, regions = regions, categories = categories, names = names
  )

  return(g)
}


#' Internal: Generate igraph Object
#'
#' @description
#' Internal helper function that generates an igraph object using specified
#' graph algorithm. Used by both simulate_igraph() and simulate_network().
#'
#' @inheritParams simulate_igraph
#'
#' @return An igraph object.
#'
#' @keywords internal
#' @noRd
.generate_graph <- function(n, model, directed, weighted, weights,
                            p, m, power, m_ba, nei, p_rewire,
                            blocks, p_within, p_between, k, radius, fw, bw,
                            name_source, regions, categories, names) {
  # Generate base graph based on model

  g <- switch(model,
    "er" = {
      if (!is.null(m)) {
        igraph::sample_gnm(n, m, directed = directed)
      } else {
        igraph::sample_gnp(n, p, directed = directed)
      }
    },
    "ba" = igraph::sample_pa(n, power = power, m = m_ba, directed = directed),
    "ws" = igraph::sample_smallworld(1, n, nei, p_rewire),
    "sbm" = {
      sizes <- rep(floor(n / blocks), blocks)
      # Handle remainder
      remainder <- n - sum(sizes)
      if (remainder > 0) {
        sizes[1:remainder] <- sizes[1:remainder] + 1
      }
      pref <- matrix(p_between, blocks, blocks)
      diag(pref) <- p_within
      igraph::sample_sbm(sum(sizes), pref, sizes, directed = directed)
    },
    "reg" = igraph::sample_k_regular(n, k, directed = directed),
    "grg" = igraph::sample_grg(n, radius),
    "ff" = igraph::sample_forestfire(n, fw, bw)
  )

  # Generate node names
  actual_n <- igraph::vcount(g)
  node_names <- if (!is.null(names)) {
    if (length(names) < actual_n) {
      stop("names must have at least ", actual_n, " elements.")
    }
    names[1:actual_n]
  } else if (name_source == "human") {
    get_global_names(actual_n, regions = regions)
  } else if (name_source == "states") {
    select_states(actual_n, primary_categories = categories)
  } else {
    paste0("V", 1:actual_n)
  }

  # Assign vertex names
  igraph::V(g)$name <- node_names

  # Add weights if requested
  if (weighted && igraph::ecount(g) > 0) {
    igraph::E(g)$weight <- runif(igraph::ecount(g), weights[1], weights[2])
  }

  # Add block attribute for SBM

  if (model == "sbm") {
    sizes <- rep(floor(n / blocks), blocks)
    remainder <- n - sum(sizes)
    if (remainder > 0) {
      sizes[1:remainder] <- sizes[1:remainder] + 1
    }
    igraph::V(g)$block <- rep(1:blocks, times = sizes)
  }

  return(g)
}
