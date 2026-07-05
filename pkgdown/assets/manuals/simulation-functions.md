# Simulation Functions Guide

This guide covers all simulation functions in Saqrlab for generating synthetic TNA data.

## Quick Reference

| Function | Output | Use Case |
|----------|--------|----------|
| `simulate_sequences()` | data.frame | Raw sequence data |
| `simulate_sequences_advanced()` | data.frame | Sequences with stable patterns |
| `simulate_tna_datasets()` | list | Complete datasets (sequences + model + probabilities) |
| `simulate_tna_networks()` | list | Multiple fitted TNA/fTNA/cTNA/aTNA models |
| `simulate_tna_network()` | tna object | Single fitted TNA model |
| `simulate_group_tna_networks()` | group_tna | Group TNA model |
| `simulate_matrix()` | matrix | Transition/frequency matrices |
| `simulate_htna()` | list | Multi-type matrices for HTNA/MLNA |
| `simulate_long_data()` | tibble | Hierarchical Actor/Group/Course data |
| `simulate_onehot_data()` | tibble | One-hot encoded sequences |
| `simulate_igraph()` | igraph | Network with graph algorithms |
| `simulate_network()` | network | statnet network object |
| `simulate_edge_list()` | data.frame | Social network edge list |

---

## Sequence Generation

### simulate_sequences()

Generate Markov chain sequences from transition probabilities.

```r
# Basic usage - auto-generates probabilities
seqs <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 6,
  seed = 42
)

# With learning state names (default)
seqs <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 6,
  use_learning_states = TRUE,
  categories = "metacognitive",
  seed = 42
)

# With custom transition matrix
trans_mat <- matrix(c(
  0.7, 0.2, 0.1,
  0.1, 0.7, 0.2,
  0.2, 0.1, 0.7
), nrow = 3, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")
init_probs <- c(A = 0.5, B = 0.3, C = 0.2)

seqs <- simulate_sequences(
  trans_matrix = trans_mat,
  init_probs = init_probs,
  n_sequences = 100,
  seq_length = 20
)

# Include generating parameters in output
result <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 5,
  include_params = TRUE,
  seed = 42
)
result$sequences      # data.frame
result$trans_matrix   # generating transition matrix
result$initial_probs  # generating initial probabilities

# With NA values
seqs <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 5,
  na_range = c(2, 5),  # 2-5 NAs per sequence
  include_na = TRUE,
  seed = 42
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_sequences` | integer | 1000 | Number of sequences to generate |
| `seq_length` | integer | 20 | Maximum length of each sequence |
| `n_states` | integer | 8 | Number of states |
| `states` | character | NULL | Custom state names |
| `use_learning_states` | logical | TRUE | Use learning verb names |
| `categories` | character | "all" | Learning categories to use |
| `trans_matrix` | matrix | NULL | Custom transition matrix |
| `init_probs` | numeric | NULL | Custom initial probabilities |
| `na_range` | integer(2) | c(0,0) | Range of NAs per sequence |
| `include_na` | logical | FALSE | Whether to include NAs |
| `include_params` | logical | FALSE | Return generating parameters |
| `seed` | integer | NULL | Random seed |

---

### simulate_sequences_advanced()

Generate sequences with stable transition patterns - useful for testing edge detection.

```r
# Define stable transitions (ground truth)
stable_transitions <- list(
  c("Plan", "Monitor"),
  c("Monitor", "Evaluate"),
  c("Evaluate", "Reflect")
)

seqs <- simulate_sequences_advanced(

  trans_matrix = trans_mat,
  init_probs = init_probs,
  seq_length = 30,
  n_sequences = 200,
  stable_transitions = stable_transitions,
  stability_prob = 0.95,        # 95% follow stable paths
  unstable_mode = "random_jump", # What happens 5% of the time
  seed = 42
)
```

**Unstable Modes:**

| Mode | Description |
|------|-------------|
| `"random_jump"` | Random transition to any state |
| `"perturb_prob"` | Add noise to transition probabilities |
| `"unlikely_jump"` | Jump to low-probability state |

---

## Dataset Generation

### simulate_tna_datasets()

Generate complete datasets with sequences, fitted models, and generating probabilities.

```r
# Generate 5 complete datasets
datasets <- simulate_tna_datasets(
  n_datasets = 5,
  n_states = 6,
  n_sequences = 100,
  seq_length = 20,
  model_type = "tna",
  use_learning_states = TRUE,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

# Access components
datasets$dataset_1$sequences        # sequence data.frame
datasets$dataset_1$model            # fitted tna model
datasets$dataset_1$transition_probs # generating transition matrix
datasets$dataset_1$initial_probs    # generating initial probabilities
datasets$dataset_1$params           # all parameters used

# Generate with advanced sequences
datasets <- simulate_tna_datasets(
  n_datasets = 3,
  n_states = 5,
  use_advanced = TRUE,
  stable_transitions = list(c("Plan", "Monitor")),
  stability_prob = 0.9,
  seed = 42
)

# Control what's included
datasets <- simulate_tna_datasets(
  n_datasets = 10,
  n_states = 4,
  include_data = TRUE,   # include sequences
  include_probs = TRUE,  # include transition_probs and initial_probs
  seed = 42
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_datasets` | integer | 1 | Number of datasets to generate |
| `n_states` | integer | 8 | Number of states |
| `n_sequences` | integer | 100 | Sequences per dataset |
| `seq_length` | integer | 30 | Maximum sequence length |
| `model_type` | character | "tna" | Model type: "tna", "ftna", "ctna", "atna" |
| `use_learning_states` | logical | TRUE | Use learning verb names |
| `categories` | character | NULL | Learning categories (random if NULL) |
| `smart_select` | logical | TRUE | Intelligent state selection |
| `use_advanced` | logical | FALSE | Use advanced sequence generation |
| `stable_transitions` | list | NULL | Stable transitions (if use_advanced) |
| `stability_prob` | numeric | 0.95 | Stability probability (if use_advanced) |
| `include_data` | logical | TRUE | Include sequences in output |
| `include_probs` | logical | TRUE | Include generating probabilities |
| `seed` | integer | NULL | Random seed |
| `verbose` | logical | TRUE | Print progress messages |

---

### simulate_tna_networks()

Generate multiple fitted TNA model objects.

```r
# 10 TNA models
tna_nets <- simulate_tna_networks(
  n_networks = 10,
  n_states = 6,
  n_sequences = 150,
  model_type = "tna",
  seed = 42
)

# Access individual networks
tna_nets$network_1$weights   # transition weights
tna_nets$network_1$initial   # initial probabilities
plot(tna_nets$network_1)     # visualize

# 5 fTNA models (filtered)
ftna_nets <- simulate_tna_networks(
  n_networks = 5,
  n_states = 8,
  model_type = "ftna",
  seed = 42
)

# 5 cTNA models (continuous-time)
ctna_nets <- simulate_tna_networks(
  n_networks = 5,
  n_states = 6,
  model_type = "ctna",
  seed = 42
)

# 5 aTNA models (absorbing)
atna_nets <- simulate_tna_networks(
  n_networks = 5,
  n_states = 6,
  model_type = "atna",
  seed = 42
)
```

**Model Types:**

| Type | Description |
|------|-------------|
| `"tna"` | Standard Transition Network Analysis |
| `"ftna"` | Filtered TNA (removes weak edges) |
| `"ctna"` | Continuous-time TNA |
| `"atna"` | Absorbing TNA |

---

### simulate_tna_network()

Generate a single fitted TNA model (convenience function).

```r
# Single TNA model
model <- simulate_tna_network(
  n_states = 6,
  n_sequences = 200,
  seq_length = 25,
  seed = 42
)

model$weights
model$initial
plot(model)
tna::centralities(model)
```

---

### simulate_group_tna_networks()

Generate group TNA models for multi-group analysis.

```r
# 5 groups, 10-15 actors each
group_model <- simulate_group_tna_networks(
  n_groups = 5,
  n_actors = c(10, 15),      # variable group sizes
  n_states = 6,
  seq_length_range = c(15, 30),
  use_learning_states = TRUE,
  categories = "group_regulation",
  seed = 42
)

# Access individual group networks
names(group_model)
group_model$Group_1$weights
```

---

## Matrix Generation

### simulate_matrix()

Generate transition, frequency, or adjacency matrices.

```r
# Transition matrix (rows sum to 1)
trans_mat <- simulate_matrix(
  n_nodes = 6,
  matrix_type = "transition",
  use_learning_states = TRUE,
  seed = 42
)

# Frequency matrix
freq_mat <- simulate_matrix(
  n_nodes = 8,
  matrix_type = "frequency",
  seed = 42
)

# Adjacency matrix (binary)
adj_mat <- simulate_matrix(
  n_nodes = 10,
  matrix_type = "adjacency",
  edge_prob = 0.3,
  seed = 42
)

# Co-occurrence matrix
cooc_mat <- simulate_matrix(
  n_nodes = 6,
  matrix_type = "co-occurrence",
  seed = 42
)
```

---

### simulate_htna() / simulate_mlna() / simulate_mtna()

Generate multi-type matrices for hierarchical/multilevel network analysis.

```r
# HTNA matrix with 3 types, 5 nodes each
htna <- simulate_htna(
  n_nodes = 5,
  n_types = 3,
  type_names = c("Metacognitive", "Cognitive", "Behavioral"),
  within_prob = 0.4,    # probability of edges within types
  between_prob = 0.15,  # probability of edges between types
  categories = c("metacognitive", "cognitive", "behavioral"),
  seed = 42
)

htna$matrix      # 15x15 transition matrix
htna$node_types  # list mapping type names to node names

# Use with tna package
# tna::plot_htna(htna$matrix, htna$node_types)
# tna::plot_mlna(htna$matrix, layers = htna$node_types)
```

---

### simulate_tna_matrix()

Convenience wrapper for HTNA/MLNA matrices with learning categories.

```r
# 5 groups x 5 nodes = 25 node matrix
net <- simulate_tna_matrix(
  nodes_per_group = 5,
  group_names = c("Metacognitive", "Cognitive", "Behavioral", "Social", "Motivational"),
  seed = 42
)

net$matrix       # 25x25 matrix
net$node_types   # node groupings
```

---

## Hierarchical Data Generation

### simulate_long_data()

Generate long-format hierarchical data with Actors, Groups, and Courses.

```r
long_data <- simulate_long_data(
  n_groups = 5,
  n_actors = c(8, 12),       # 8-12 actors per group
  n_courses = 3,
  n_states = 6,
  seq_length_range = c(10, 25),
  use_learning_states = TRUE,
  seed = 42
)

head(long_data)
# Actor  Group   Course  Time  Action
# A1_G1  Group_1 Course_1  1   Plan
# A1_G1  Group_1 Course_1  2   Monitor
# ...
```

---

### simulate_onehot_data()

Generate one-hot encoded sequence data.

```r
onehot_data <- simulate_onehot_data(
  n_groups = 3,
  n_actors = 10,
  n_states = 5,
  seed = 42
)

head(onehot_data)
# Actor  Group  Time  Plan  Monitor  Evaluate  Reflect  Regulate
# A1_G1  G1     1     1     0        0         0        0
# A1_G1  G1     2     0     1        0         0        0
# ...
```

---

## Network Object Generation

### simulate_igraph()

Generate igraph network objects with various graph models.

```r
# Barabasi-Albert preferential attachment
g <- simulate_igraph(
  n = 30,
  model = "ba",
  directed = TRUE,
  weighted = TRUE,
  name_source = "human",
  regions = "arab",
  seed = 42
)

# Watts-Strogatz small world
g <- simulate_igraph(
  n = 50,
  model = "ws",
  k = 4,        # each node connected to k neighbors
  p = 0.1,      # rewiring probability
  seed = 42
)

# Stochastic block model
g <- simulate_igraph(
  n = 40,
  model = "sbm",
  blocks = 4,
  seed = 42
)
```

**Graph Models:**

| Model | Description |
|-------|-------------|
| `"er"` | Erdos-Renyi random graph |
| `"ba"` | Barabasi-Albert preferential attachment |
| `"ws"` | Watts-Strogatz small world |
| `"sbm"` | Stochastic block model |
| `"reg"` | Regular graph (fixed degree) |
| `"grg"` | Geometric random graph |
| `"ff"` | Forest fire model |

---

### simulate_network()

Generate statnet network objects.

```r
net <- simulate_network(
  n = 30,
  model = "ba",
  directed = TRUE,
  name_source = "human",
  regions = "europe",
  seed = 42
)
```

---

### simulate_edge_list()

Generate social network edge lists.

```r
edges <- simulate_edge_list(
  n_nodes = 20,
  n_edges = 50,
  directed = TRUE,
  weighted = TRUE,
  name_source = "human",
  regions = "africa",
  seed = 42
)

head(edges)
# from    to      weight
# Amara   Kwame   0.73
# ...
```

---

## Utility Functions

### generate_probabilities()

Generate random transition matrices and initial probabilities.

```r
probs <- generate_probabilities(
  n_states = 5,
  seed = 42
)

probs$transition_matrix  # 5x5 matrix (rows sum to 1)
probs$initial_probs      # 5 probabilities (sum to 1)
```

---

### generate_param_grid()

Create parameter grids for simulation studies.

```r
# Latin Hypercube Sampling
grid <- generate_param_grid(
  param_ranges = list(
    n_sequences = c(50, 500),
    seq_length = c(15, 50),
    n_states = c(4, 10)
  ),
  n = 100,
  method = "lhs"
)

# Regular grid
grid <- generate_param_grid(
  param_ranges = list(
    n_sequences = c(100, 200, 300),
    seq_length = c(20, 30, 40)
  ),
  n = 20,
  method = "grid"
)
```

---

### select_states()

Intelligently select learning states based on network size.

```r
# Small network: single category
states_5 <- select_states(5, seed = 42)

# Medium network: balanced categories
states_10 <- select_states(10, seed = 42)

# Large network: multiple categories
states_20 <- select_states(20, seed = 42)

# Biased selection
srl_states <- select_states(
  n_states = 10,
  primary_categories = "metacognitive",
  secondary_categories = "cognitive",
  primary_ratio = 0.7,
  seed = 42
)
```

---

## Learning States

Saqrlab includes 180+ learning action verbs organized into 8 categories:

| Category | Examples | Count |
|----------|----------|-------|
| `metacognitive` | Plan, Monitor, Evaluate, Reflect | 20 |
| `cognitive` | Read, Analyze, Summarize, Memorize | 30 |
| `behavioral` | Practice, Annotate, Research, Write | 30 |
| `social` | Collaborate, Discuss, Explain, Share | 30 |
| `motivational` | Focus, Persist, Explore, Strive | 30 |
| `affective` | Enjoy, Appreciate, Cope, Manage | 30 |
| `group_regulation` | Adapt, Cohesion, Consensus, Coregulate | 9 |
| `lms` | View, Access, Download, Submit, Quiz | 30 |

```r
# Get all metacognitive verbs
get_learning_states("metacognitive")

# Get 10 random verbs from cognitive + behavioral
get_learning_states(c("cognitive", "behavioral"), n = 10, seed = 42)

# View all categories
list_learning_categories()
```

---

## Global Names

300 diverse human names from multiple regions:

```r
# Arab names
get_global_names(10, regions = "arab")

# European names
get_global_names(10, regions = "europe")

# Mixed regions
get_global_names(20, regions = c("arab", "east_asia", "africa"))

# View all regions
list_name_regions()
```

**Region Shortcuts:** `europe`, `middle_east`, `asia`, `africa`, `americas`, `oceania`

---

## Reproducibility

All simulation functions accept a `seed` parameter:

```r
# Same seed = identical results
mat1 <- simulate_matrix(n_nodes = 5, seed = 123)
mat2 <- simulate_matrix(n_nodes = 5, seed = 123)
identical(mat1, mat2)  # TRUE

# Different seeds = different results
mat3 <- simulate_matrix(n_nodes = 5, seed = 456)
identical(mat1, mat3)  # FALSE
```

---

## Backward Compatibility

All old function names still work (silent aliases):

| Old Name | New Name |
|----------|----------|
| `generate_tna_datasets()` | `simulate_tna_datasets()` |
| `generate_tna_networks()` | `simulate_tna_networks()` |
| `generate_group_tna_networks()` | `simulate_group_tna_networks()` |
| `generate_tna_matrix()` | `simulate_tna_matrix()` |
| `generate_sequence_data()` | `simulate_tna_datasets()` |
| `analyze_grid_results()` | `summarize_grid_results()` |
| `calculate_edge_recovery()` | `compare_edge_recovery()` |
| `evaluate_bootstrap()` | `run_bootstrap_iteration()` |
| `smart_select_states()` | `select_states()` |
| `create_param_grid()` | `generate_param_grid()` |
