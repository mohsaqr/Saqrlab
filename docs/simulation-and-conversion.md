# Simulation and Conversion Functions

This tutorial covers all functions in Saqrlab that generate data or convert between formats.

## Quick Reference

### Simulation Functions

| Function | Output | Description |
|----------|--------|-------------|
| `simulate_igraph()` | igraph | Network object with graph algorithms |
| `simulate_network()` | network | statnet network object |
| `simulate_tna_network()` | tna model | Ready-to-use fitted TNA model |
| `simulate_matrix()` | matrix | Transition/frequency/adjacency matrix |
| `simulate_htna()` | list | Multi-type matrix for HTNA/MLNA |
| `simulate_sequences()` | data.frame | Markov chain sequences |
| `simulate_sequences_advanced()` | data.frame | Sequences with stable patterns |
| `simulate_long_data()` | tibble | Hierarchical Actor/Group/Course data |
| `simulate_onehot_data()` | tibble | One-hot encoded sequences |
| `simulate_edge_list()` | data.frame | Social network edge list |

### Generation Functions

| Function | Output | Description |
|----------|--------|-------------|
| `generate_probabilities()` | list | Random transition + initial probabilities |
| `generate_tna_datasets()` | list | Complete datasets with sequences + model |
| `generate_tna_networks()` | list | Multiple fitted TNA models |
| `generate_group_tna_networks()` | group_tna | Group TNA model |
| `generate_tna_matrix()` | list | Matrix with node groupings |

### Conversion Functions

| Function | Input | Output | Description |
|----------|-------|--------|-------------|
| `wide_to_long()` | wide df | long df | One row per sequence → one row per action |
| `long_to_wide()` | long df | wide df | One row per action → one row per sequence |
| `prepare_for_tna()` | any df | wide df | Prepare data for tna::tna() |
| `action_to_onehot()` | long df | onehot df | Action column → binary columns |

### Extraction Functions

| Function | Input | Output | Description |
|----------|-------|--------|-------------|
| `extract_transition_matrix()` | model | matrix | Get transition matrix |
| `extract_initial_probs()` | model | vector | Get initial probabilities |
| `extract_edges()` | model | data.frame | Get edge list |

---

## 1. Network Objects

### simulate_igraph() - igraph Networks

Generate network objects using common graph algorithms with human or learning state names.

```r
library(Saqrlab)
library(igraph)

# Simplest: random network with 20-50 nodes
g <- simulate_igraph(seed = 42)
vcount(g)  # Random between 20-50
V(g)$name  # Human names from diverse regions

# Specify size
g <- simulate_igraph(n = 30, seed = 42)
```

**Graph Algorithms:**

```r
# Erdos-Renyi (random) - default
g_er <- simulate_igraph(n = 30, model = "er", p = 0.1, seed = 42)

# Barabasi-Albert (scale-free)
g_ba <- simulate_igraph(n = 30, model = "ba", m_ba = 2, seed = 42)

# Watts-Strogatz (small-world)
g_ws <- simulate_igraph(n = 30, model = "ws", nei = 3, p_rewire = 0.1, seed = 42)

# Stochastic Block Model (communities)
g_sbm <- simulate_igraph(n = 30, model = "sbm", blocks = 3, seed = 42)
V(g_sbm)$block  # Community assignments

# Regular graph (fixed degree)
g_reg <- simulate_igraph(n = 30, model = "reg", k = 4, seed = 42)

# Geometric random graph (spatial)
g_grg <- simulate_igraph(n = 30, model = "grg", radius = 0.25, seed = 42)

# Forest fire (growing network)
g_ff <- simulate_igraph(n = 30, model = "ff", seed = 42)
```

**Name Sources:**

```r
# Human names from specific regions
g <- simulate_igraph(n = 20, regions = "arab", seed = 42)
g <- simulate_igraph(n = 20, regions = "africa", seed = 42)
g <- simulate_igraph(n = 20, regions = c("east_asia", "south_asia"), seed = 42)

# Learning states instead of human names
g <- simulate_igraph(n = 15, name_source = "states", seed = 42)
V(g)$name  # "Plan", "Monitor", "Evaluate", etc.

# Custom names
g <- simulate_igraph(n = 5, names = c("A", "B", "C", "D", "E"), seed = 42)
```

**Weighted Networks:**

```r
g <- simulate_igraph(n = 20, weighted = TRUE, weights = c(0.1, 1.0), seed = 42)
E(g)$weight
```

### simulate_network() - statnet Networks

Same as simulate_igraph() but returns statnet network objects.

```r
library(network)
library(sna)

net <- simulate_network(n = 30, model = "ba", seed = 42)
network.size(net)
network.vertex.names(net)

# Use with sna functions
betweenness(net)
closeness(net)
```

### simulate_tna_network() - Fitted TNA Models

The simplest way to get a ready-to-use TNA model.

```r
# Generate a fitted tna model
model <- simulate_tna_network(seed = 42)
class(model)  # "tna"

# Use with tna package
library(tna)
plot(model)
centralities(model)
communities(model)

# Custom configuration
model <- simulate_tna_network(
  n_states = 8,
  n_sequences = 500,
  categories = "group_regulation",
  seed = 123
)
```

---

## 2. Matrices

### simulate_matrix() - Transition Matrices

Generate transition matrices with learning state names.

```r
# Default: 9-node transition matrix
mat <- simulate_matrix(seed = 42)
dim(mat)
rowSums(mat)  # All rows sum to 1
rownames(mat)  # Learning state names
```

**Matrix Types:**

```r
# Transition (default) - rows sum to 1
trans_mat <- simulate_matrix(n_nodes = 6, matrix_type = "transition", seed = 42)

# Frequency - integer counts
freq_mat <- simulate_matrix(n_nodes = 6, matrix_type = "frequency", seed = 42)

# Co-occurrence - symmetric
cooc_mat <- simulate_matrix(n_nodes = 6, matrix_type = "co-occurrence", seed = 42)
isSymmetric(cooc_mat)  # TRUE

# Adjacency - binary or weighted
adj_mat <- simulate_matrix(n_nodes = 6, matrix_type = "adjacency", seed = 42)
```

### simulate_htna() - Multi-Type Matrices

Generate matrices with node type groupings for HTNA/MLNA analysis.

```r
# Default: 5 types x 5 nodes = 25-node matrix
net <- simulate_htna(seed = 42)
dim(net$matrix)           # 25 x 25
names(net$node_types)     # Metacognitive, Cognitive, Behavioral, Social, Motivational
net$node_types$Metacognitive  # Nodes in this type

# Custom configuration
net <- simulate_htna(
  n_nodes = 4,           # Nodes per type
  n_types = 3,           # Number of types
  type_names = c("Macro", "Meso", "Micro"),
  within_prob = 0.5,     # Higher within-type connectivity
  between_prob = 0.1,    # Lower between-type connectivity
  seed = 42
)

# Use with tna package
library(tna)
plot_htna(net$matrix, net$node_types, layout = "polygon")
plot_mlna(net$matrix, layers = net$node_types)
```

Aliases: `simulate_mlna()`, `simulate_mtna()` (identical functions)

---

## 3. Sequences

### simulate_sequences() - Basic Sequences

Generate Markov chain sequences.

```r
# Auto-generate with learning states
sequences <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 6,
  seed = 42
)
head(sequences)
dim(sequences)  # 100 rows x 20 columns
```

**With Custom Probabilities:**

```r
# Define your own transition matrix
trans_mat <- matrix(c(
  0.7, 0.2, 0.1,
  0.2, 0.5, 0.3,
  0.1, 0.3, 0.6
), nrow = 3, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("Plan", "Execute", "Review")

init_probs <- c(Plan = 0.5, Execute = 0.3, Review = 0.2)

sequences <- simulate_sequences(
  trans_matrix = trans_mat,
  init_probs = init_probs,
  n_sequences = 100,
  seq_length = 15
)
```

**Learning Categories:**

```r
# Use specific categories
sequences <- simulate_sequences(
  n_sequences = 50,
  seq_length = 20,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)
unique(unlist(sequences))  # States from those categories
```

**Variable Lengths (NAs):**

```r
# Add 0-5 trailing NAs per sequence
sequences <- simulate_sequences(
  n_sequences = 50,
  seq_length = 20,
  n_states = 5,
  na_range = c(0, 5),
  include_na = TRUE,
  seed = 42
)
```

**Get Generating Parameters:**

```r
result <- simulate_sequences(
  n_sequences = 50,
  seq_length = 15,
  n_states = 4,
  include_params = TRUE,
  seed = 42
)
result$sequences
result$trans_matrix
result$initial_probs
result$state_names
```

### simulate_sequences_advanced() - Stable Patterns

Generate sequences with stable transition patterns.

```r
# Define stable transitions
stable_transitions <- list(
  c("Plan", "Monitor"),
  c("Monitor", "Evaluate")
)

# Create matrix with these states
trans_mat <- matrix(runif(16), 4, 4)
trans_mat <- trans_mat / rowSums(trans_mat)
rownames(trans_mat) <- colnames(trans_mat) <- c("Plan", "Monitor", "Evaluate", "Execute")
init_probs <- c(Plan = 0.4, Monitor = 0.2, Evaluate = 0.2, Execute = 0.2)

# Generate with 90% stability
sequences <- simulate_sequences_advanced(
  trans_matrix = trans_mat,
  init_probs = init_probs,
  n_sequences = 100,
  seq_length = 30,
  stable_transitions = stable_transitions,
  stability_prob = 0.90,
  seed = 42
)
```

---

## 4. Hierarchical Data

### simulate_long_data() - Long Format

Generate hierarchical data with Actor/Group/Course structure.

```r
# Educational data: 5 groups, 10 actors each, 3 courses
long_data <- simulate_long_data(
  n_groups = 5,
  n_actors = 10,
  n_courses = 3,
  categories = "group_regulation",
  seq_length_range = c(10, 25),
  seed = 42
)
head(long_data)
# Columns: Actor, Achiever, Group, Course, Time, Action
```

**Variable Group Sizes:**

```r
long_data <- simulate_long_data(
  n_groups = 5,
  n_actors = c(8, 12),  # Min and max
  n_courses = 2,
  seed = 42
)
```

**Achievement Levels:**

```r
long_data <- simulate_long_data(
  n_groups = 5,
  n_actors = 10,
  achiever_levels = c("High", "Medium", "Low"),
  achiever_probs = c(0.3, 0.5, 0.2),
  seed = 42
)
table(long_data$Achiever)
```

### simulate_onehot_data() - One-Hot Encoding

Generate one-hot encoded hierarchical data.

```r
onehot_data <- simulate_onehot_data(
  n_groups = 3,
  n_actors = 10,
  n_states = 4,
  seed = 42
)
head(onehot_data)
# Columns: Actor, Group, Time, State1, State2, State3, State4 (0/1 values)
```

### simulate_edge_list() - Social Networks

Generate edge lists for social network analysis.

```r
edges <- simulate_edge_list(
  n_nodes = 20,
  n_edges = 50,
  directed = TRUE,
  n_classes = 3,
  seed = 42
)
head(edges)
# Columns: source, target, weight, class
```

---

## 5. Generation Functions

### generate_probabilities() - Random Probabilities

Generate transition and initial probabilities.

```r
probs <- generate_probabilities(n_states = 5, seed = 42)
probs$initial_probs      # Named vector summing to 1
probs$transition_probs   # Matrix with rows summing to 1
probs$state_names        # A, B, C, D, E

# With custom state names
probs <- generate_probabilities(
  n_states = 4,
  states = c("Plan", "Execute", "Review", "Complete"),
  seed = 42
)
```

### generate_tna_datasets() - Complete Datasets

Generate datasets with sequences, fitted models, and generating probabilities.

```r
datasets <- generate_tna_datasets(
  n_datasets = 5,
  n_states = 6,
  n_sequences = 100,
  seq_length = 20,
  seed = 42
)

# Access components
datasets$dataset_1$sequences       # Sequence data
datasets$dataset_1$model           # Fitted TNA model
datasets$dataset_1$transition_probs # Generating transition matrix
datasets$dataset_1$initial_probs   # Generating initial probabilities
datasets$dataset_1$params          # All parameters
```

### generate_tna_networks() - Multiple TNA Models

Generate multiple fitted TNA models.

```r
networks <- generate_tna_networks(
  n_networks = 5,
  n_states = 6,
  n_sequences = 150,
  model_type = "tna",
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

# Use the networks
library(tna)
plot(networks$network_1)
centralities(networks$network_1)
```

### generate_group_tna_networks() - Group TNA

Generate group TNA models with multiple groups.

```r
group_model <- generate_group_tna_networks(
  n_groups = 4,
  n_actors = 15,
  n_states = 5,
  seq_length_range = c(10, 30),
  seed = 42
)

# Access individual group networks
names(group_model)
```

### generate_tna_matrix() - Matrix with Groupings

Generate matrix with node groupings for HTNA/MLNA plots.

```r
net <- generate_tna_matrix(
  nodes_per_group = 5,
  group_names = c("Macro", "Meso", "Micro"),
  seed = 42
)
net$matrix       # Transition matrix
net$node_types   # Node groupings for plotting
```

---

## 6. Conversion Functions

### wide_to_long() - Wide to Long Format

Convert sequence data from wide format (rows = sequences) to long format (rows = actions).

```r
# Start with wide format
sequences <- simulate_sequences(n_sequences = 50, seq_length = 10, seed = 42)
head(sequences)
# V1, V2, V3, ..., V10

# Convert to long format
long_data <- wide_to_long(sequences)
head(long_data)
# id, Time, Action

# With custom column names
long_data <- wide_to_long(
  sequences,
  id_col = NULL,           # Auto-generate IDs
  time_prefix = "V",       # Columns to pivot
  action_col = "State",    # Output column name
  time_col = "Step",       # Output column name
  drop_na = TRUE           # Remove NA values
)
```

### long_to_wide() - Long to Wide Format

Convert long format (rows = actions) to wide format (rows = sequences).

```r
# Start with long format
long_data <- simulate_long_data(n_groups = 3, n_actors = 10, seed = 42)
head(long_data)

# Convert to wide format
wide_data <- long_to_wide(
  long_data,
  id_col = "Actor",
  time_col = "Time",
  action_col = "Action"
)
head(wide_data)
# Actor, V1, V2, V3, ...
```

### prepare_for_tna() - Prepare for TNA

Prepare any data format for use with tna::tna().

```r
# From wide format
sequences <- simulate_sequences(n_sequences = 100, seed = 42)
tna_data <- prepare_for_tna(sequences, type = "sequences")
model <- tna::tna(tna_data)

# From long format
long_data <- simulate_long_data(n_groups = 5, seed = 42)
tna_data <- prepare_for_tna(long_data, type = "long", id_col = "Actor")
model <- tna::tna(tna_data)

# Auto-detect format
tna_data <- prepare_for_tna(some_data, type = "auto")
```

### action_to_onehot() - One-Hot Encoding

Convert categorical Action column to binary indicator columns.

```r
# Start with long format
long_data <- simulate_long_data(n_groups = 3, n_states = 4, seed = 42)
head(long_data)

# Convert to one-hot
onehot_data <- action_to_onehot(long_data)
head(onehot_data)
# Actor, Group, Time, Plan, Monitor, Execute, Review (0/1 values)

# With prefix
onehot_data <- action_to_onehot(long_data, prefix = "state_")
# state_Plan, state_Monitor, etc.
```

---

## 7. Extraction Functions

### extract_transition_matrix() - Get Transition Matrix

Extract transition matrix from a TNA model.

```r
model <- simulate_tna_network(seed = 42)

# Raw weights
trans_mat <- extract_transition_matrix(model, type = "raw")

# Row-normalized (rows sum to 1)
trans_mat <- extract_transition_matrix(model, type = "scaled")
rowSums(trans_mat)  # All 1
```

### extract_initial_probs() - Get Initial Probabilities

Extract initial state probabilities from a TNA model.

```r
model <- simulate_tna_network(seed = 42)
init_probs <- extract_initial_probs(model)
sum(init_probs)  # 1
```

### extract_edges() - Get Edge List

Extract edge list from a TNA model.

```r
model <- simulate_tna_network(seed = 42)

edges <- extract_edges(
  model,
  threshold = 0.05,      # Minimum weight to include
  include_self = FALSE,  # Exclude self-loops
  sort_by = "weight"     # Sort by weight descending
)
head(edges)
# from, to, weight

# Use with igraph
library(igraph)
g <- graph_from_data_frame(edges, directed = TRUE)
plot(g)
```

---

## 8. Workflow Examples

### Example 1: Generate and Analyze a Network

```r
library(Saqrlab)
library(igraph)

# Generate network
g <- simulate_igraph(n = 30, model = "ba", regions = "europe", seed = 42)

# Analyze
degree(g)
betweenness(g)
closeness(g)

# Plot
plot(g, vertex.size = 10, vertex.label.cex = 0.7)
```

### Example 2: Simulate TNA Workflow

```r
library(Saqrlab)
library(tna)

# Generate sequences
sequences <- simulate_sequences(
  n_sequences = 200,
  seq_length = 25,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

# Fit TNA model
model <- tna(sequences)

# Analyze
plot(model)
centralities(model)
communities(model)
```

### Example 3: Convert Between Formats

```r
library(Saqrlab)

# Generate long format data
long_data <- simulate_long_data(
  n_groups = 5,
  n_actors = 15,
  categories = "group_regulation",
  seed = 42
)

# Convert to wide for TNA
wide_data <- long_to_wide(long_data, id_col = "Actor")

# Prepare for TNA
tna_data <- prepare_for_tna(wide_data, type = "sequences")

# Fit model
library(tna)
model <- tna(tna_data)
```

### Example 4: Extract and Recreate

```r
library(Saqrlab)

# Fit a model
model <- simulate_tna_network(n_states = 5, seed = 42)

# Extract components
trans_mat <- extract_transition_matrix(model, type = "scaled")
init_probs <- extract_initial_probs(model)

# Generate new sequences from the same parameters
new_sequences <- simulate_sequences(
  trans_matrix = trans_mat,
  init_probs = init_probs,
  n_sequences = 500,
  seq_length = 30
)
```

---

## 9. Name Sources

### Human Names (GLOBAL_NAMES)

1,472 culturally diverse names from 25 regions.

```r
# See all regions
list_name_regions()

# Get names from specific regions
get_global_names(20, regions = "arab", seed = 42)
get_global_names(20, regions = "east_asia", seed = 42)

# Use shortcuts
get_global_names(20, regions = "europe", seed = 42)     # All European regions
get_global_names(20, regions = "africa", seed = 42)      # All African regions
get_global_names(20, regions = "asia", seed = 42)        # All Asian regions
```

### Learning States (LEARNING_STATES)

180+ learning action verbs in 8 categories.

```r
# See all categories
list_learning_categories()

# Get states from specific categories
get_learning_states("metacognitive")
get_learning_states("cognitive", n = 10)
get_learning_states(c("metacognitive", "cognitive"), n = 15)

# Smart selection based on n_states
smart_select_states(n_states = 10, seed = 42)
```
