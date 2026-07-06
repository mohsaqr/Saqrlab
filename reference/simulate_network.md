# Simulate statnet Network Object

Generate a statnet network object using common graph algorithms with
realistic node names from human names or learning states.

## Usage

``` r
simulate_network(
  n = NULL,
  model = c("er", "ba", "ws", "sbm", "reg", "grg", "ff"),
  name_source = c("human", "states"),
  regions = "all",
  categories = "all",
  names = NULL,
  directed = FALSE,
  weighted = FALSE,
  weights = c(0.1, 1),
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
  seed = NULL
)
```

## Arguments

- n:

  Integer or NULL. Number of nodes. If NULL (default), randomly selects
  between 20-50 nodes.

- model:

  Character. Graph generation algorithm:

  - `"er"`: Erdos-Renyi random graph

  - `"ba"`: Barabasi-Albert scale-free network

  - `"ws"`: Watts-Strogatz small-world network

  - `"sbm"`: Stochastic Block Model (community structure)

  - `"reg"`: Regular graph (fixed degree)

  - `"grg"`: Geometric Random Graph (spatial)

  - `"ff"`: Forest Fire (growing network)

  Default: "er".

- name_source:

  Character. Source for node names:

  - `"human"`: Culturally diverse human names from GLOBAL_NAMES

  - `"states"`: Learning action verbs from LEARNING_STATES

  Default: "human".

- regions:

  Character vector. Regions to sample human names from (only used when
  name_source = "human"). Can be specific regions (e.g., "arab",
  "east_asia"), shortcuts (e.g., "europe", "africa", "asia"), or "all".
  See
  [`list_name_regions`](https://pak.dynasite.org/Saqrlab/reference/list_name_regions.md).
  Default: "all".

- categories:

  Character vector. Learning state categories (only used when
  name_source = "states"). Options: "metacognitive", "cognitive",
  "behavioral", "social", "motivational", "affective",
  "group_regulation", "lms", or "all". Default: "all".

- names:

  Character vector or NULL. Custom node names. Overrides name_source if
  provided. Default: NULL.

- directed:

  Logical. If TRUE, generate directed network. Default: FALSE.

- weighted:

  Logical. If TRUE, add random edge weights. Default: FALSE.

- weights:

  Numeric vector of length 2. Weight range \[min, max\]. Default: c(0.1,
  1.0).

- p:

  Numeric. Edge probability for Erdos-Renyi model. Default: 0.1.

- m:

  Integer or NULL. Fixed number of edges for Erdos-Renyi. Overrides p if
  provided. Default: NULL.

- power:

  Numeric. Attachment power for Barabasi-Albert. Default: 1.

- m_ba:

  Integer. Edges per new vertex for Barabasi-Albert. Default: 2.

- nei:

  Integer. Neighborhood size for Watts-Strogatz. Default: 2.

- p_rewire:

  Numeric. Rewiring probability for Watts-Strogatz. Default: 0.05.

- blocks:

  Integer. Number of blocks for SBM. Default: 3.

- p_within:

  Numeric. Within-block edge probability for SBM. Default: 0.3.

- p_between:

  Numeric. Between-block edge probability for SBM. Default: 0.05.

- k:

  Integer. Degree for regular graphs. Default: 4.

- radius:

  Numeric. Connection radius for geometric random graph. Default: 0.25.

- fw:

  Numeric. Forward burning probability for Forest Fire. Default: 0.35.

- bw:

  Numeric. Backward burning factor for Forest Fire. Default: 0.32.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

## Value

A network object (class "network") with vertex names and optional edge
weights.

## Details

This function generates networks using igraph algorithms internally,
then converts to statnet's network class. The resulting object is
compatible with all sna and network package functions.

Vertex names are stored in the "vertex.names" attribute and can be
accessed with
[`network.vertex.names()`](https://rdrr.io/pkg/network/man/attribute.methods.html).

## See also

[`simulate_igraph`](https://pak.dynasite.org/Saqrlab/reference/simulate_igraph.md)
for igraph objects,
[`simulate_matrix`](https://pak.dynasite.org/Saqrlab/reference/simulate_matrix.md)
for adjacency matrices,
[`simulate_edge_list`](https://pak.dynasite.org/Saqrlab/reference/simulate_edge_list.md)
for edge list data frames.

## Examples

``` r
if (FALSE) { # \dontrun{
library(network)
library(sna)

# Default: Erdos-Renyi with human names
net <- simulate_network(n = 20, seed = 42)
class(net)  # "network"
network.vertex.names(net)  # Diverse human names

# Names from specific regions
net_asia <- simulate_network(n = 20, regions = "asia", seed = 42)
network.vertex.names(net_asia)  # Asian names

# Use learning state names instead
net_states <- simulate_network(n = 15, name_source = "states", seed = 42)
network.vertex.names(net_states)  # Action verbs

# Scale-free for SNA analysis
net_sf <- simulate_network(n = 50, model = "ba", seed = 42)
betweenness(net_sf)
closeness(net_sf)

# Community structure
net_sbm <- simulate_network(n = 30, model = "sbm", blocks = 3, seed = 42)

# Weighted network
net_w <- simulate_network(n = 20, model = "er", weighted = TRUE, seed = 42)
net_w %e% "weight"  # Edge weights
} # }
```
