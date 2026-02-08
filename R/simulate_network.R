#' Simulate statnet Network Object
#'
#' @description
#' Generate a statnet network object using common graph algorithms with
#' realistic node names from human names or learning states.
#'
#' @inheritParams simulate_igraph
#'
#' @return A network object (class "network") with vertex names and optional
#'   edge weights.
#'
#' @details
#' This function generates networks using igraph algorithms internally, then
#' converts to statnet's network class. The resulting object is compatible
#' with all sna and network package functions.
#'
#' Vertex names are stored in the "vertex.names" attribute and can be accessed
#' with \code{network.vertex.names()}.
#'
#' @examples
#' \dontrun{
#' library(network)
#' library(sna)
#'
#' # Default: Erdos-Renyi with human names
#' net <- simulate_network(n = 20, seed = 42)
#' class(net)  # "network"
#' network.vertex.names(net)  # Diverse human names
#'
#' # Names from specific regions
#' net_asia <- simulate_network(n = 20, regions = "asia", seed = 42)
#' network.vertex.names(net_asia)  # Asian names
#'
#' # Use learning state names instead
#' net_states <- simulate_network(n = 15, name_source = "states", seed = 42)
#' network.vertex.names(net_states)  # Action verbs
#'
#' # Scale-free for SNA analysis
#' net_sf <- simulate_network(n = 50, model = "ba", seed = 42)
#' betweenness(net_sf)
#' closeness(net_sf)
#'
#' # Community structure
#' net_sbm <- simulate_network(n = 30, model = "sbm", blocks = 3, seed = 42)
#'
#' # Weighted network
#' net_w <- simulate_network(n = 20, model = "er", weighted = TRUE, seed = 42)
#' net_w %e% "weight"  # Edge weights
#' }
#'
#' @seealso
#' \code{\link{simulate_igraph}} for igraph objects,
#' \code{\link{simulate_matrix}} for adjacency matrices,
#' \code{\link{simulate_edge_list}} for edge list data frames.
#'
#' @export
simulate_network <- function(n = NULL,
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
  # Check if network package is available
  if (!requireNamespace("network", quietly = TRUE)) {
    stop("Package 'network' is required for simulate_network(). ",
         "Install it with: install.packages('network')")
  }

  # Match arguments
  model <- match.arg(model)
  name_source <- match.arg(name_source)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Random n if not specified
  if (is.null(n)) {
    n <- sample(20:50, 1)
  }

  # Generate igraph object using internal helper
  g <- .generate_graph(
    n = n, model = model, directed = directed, weighted = weighted,
    weights = weights, p = p, m = m, power = power, m_ba = m_ba,
    nei = nei, p_rewire = p_rewire, blocks = blocks, p_within = p_within,
    p_between = p_between, k = k, radius = radius, fw = fw, bw = bw,
    name_source = name_source, regions = regions, categories = categories, names = names
  )

  # Convert to network object
  adj <- igraph::as_adjacency_matrix(g, sparse = FALSE, attr = if (weighted) "weight" else NULL)
  net <- network::network(adj, directed = igraph::is_directed(g),
                          ignore.eval = !weighted,
                          names.eval = if (weighted) "weight" else NULL)

  # Transfer vertex names
  network::set.vertex.attribute(net, "vertex.names", igraph::V(g)$name)

  # Transfer block attribute for SBM
  if (model == "sbm") {
    network::set.vertex.attribute(net, "block", igraph::V(g)$block)
  }

  return(net)
}
