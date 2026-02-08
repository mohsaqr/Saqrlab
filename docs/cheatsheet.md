# Saqrlab Cheatsheet

## Network Objects

```r
# igraph network (20-50 random nodes by default)
g <- simulate_igraph(seed = 42)
g <- simulate_igraph(n = 30, model = "ba", regions = "arab", seed = 42)

# statnet network
net <- simulate_network(n = 30, model = "sbm", blocks = 3, seed = 42)

# Fitted TNA model
model <- simulate_tna_network(n_states = 8, seed = 42)
```

**Graph Models:** `"er"` (Erdos-Renyi), `"ba"` (Barabasi-Albert), `"ws"` (Watts-Strogatz), `"sbm"` (Stochastic Block), `"reg"` (Regular), `"grg"` (Geometric), `"ff"` (Forest Fire)

**Name Sources:** `name_source = "human"` (default) or `"states"` (learning verbs)

**Regions:** `"arab"`, `"east_asia"`, `"nordic"`, `"africa"`, `"europe"`, `"asia"`, etc.

---

## Matrices

```r
# Transition matrix (rows sum to 1)
mat <- simulate_matrix(n_nodes = 9, seed = 42)

# Multi-type matrix for HTNA/MLNA
net <- simulate_htna(n_nodes = 5, n_types = 3, seed = 42)
# net$matrix, net$node_types
```

**Matrix Types:** `"transition"`, `"frequency"`, `"co-occurrence"`, `"adjacency"`

---

## Sequences

```r
# Basic sequences
seq <- simulate_sequences(n_sequences = 100, seq_length = 20, n_states = 6, seed = 42)

# With custom probabilities
seq <- simulate_sequences(trans_matrix = mat, init_probs = probs, n_sequences = 100)

# Get parameters back
result <- simulate_sequences(..., include_params = TRUE)
# result$sequences, result$trans_matrix, result$initial_probs
```

---

## Hierarchical Data

```r
# Long format (Actor/Group/Course)
long <- simulate_long_data(n_groups = 5, n_actors = 10, n_courses = 3, seed = 42)

# One-hot encoded
onehot <- simulate_onehot_data(n_groups = 3, n_actors = 10, seed = 42)

# Edge list
edges <- simulate_edge_list(n_nodes = 20, n_edges = 50, seed = 42)
```

---

## Generation

```r
# Random probabilities
probs <- generate_probabilities(n_states = 5, seed = 42)

# Complete datasets (sequences + model + probabilities)
data <- simulate_tna_datasets(n_datasets = 5, n_states = 6, seed = 42)

# Multiple TNA models
nets <- simulate_tna_networks(n_networks = 5, n_states = 6, seed = 42)

# Group TNA
group_net <- simulate_group_tna_networks(n_groups = 4, n_actors = 15, seed = 42)
```

---

## Conversion

```r
# Wide to Long
long <- wide_to_long(sequences)

# Long to Wide
wide <- long_to_wide(long_data, id_col = "Actor")

# Prepare for TNA
tna_data <- prepare_for_tna(data, type = "auto")

# One-hot encoding
onehot <- action_to_onehot(long_data)
```

---

## Extraction

```r
# From TNA model
trans_mat <- extract_transition_matrix(model, type = "scaled")
init_probs <- extract_initial_probs(model)
edges <- extract_edges(model, threshold = 0.05)
```

---

## Names

```r
# Human names by region
get_global_names(20, regions = "arab")
get_global_names(20, regions = "africa")  # Shortcut for all African regions
list_name_regions()  # See all regions

# Learning states by category
get_learning_states("metacognitive")
get_learning_states(c("metacognitive", "cognitive"), n = 10)
select_states(n_states = 10)  # Auto-select based on size
list_learning_categories()  # See all categories
```

**Region Shortcuts:** `europe`, `middle_east`, `asia`, `africa`, `americas`, `oceania`

**Learning Categories:** `metacognitive`, `cognitive`, `behavioral`, `social`, `motivational`, `affective`, `group_regulation`, `lms`
