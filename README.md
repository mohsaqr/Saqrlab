# Saqrlab <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/mohsaqr/Saqrlab/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mohsaqr/Saqrlab/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Simulation and Analysis Tools for Temporal Network Analysis (TNA)**

Saqrlab is an R package providing comprehensive tools for simulating, analyzing, and comparing Temporal Network Analysis models. Designed for educational researchers, it includes:

- Markov chain sequence simulation with realistic learning state names
- Multiple TNA model fitting (TNA, fTNA, cTNA, aTNA)
- Network comparison metrics (correlation, RMSE, edge recovery)
- Bootstrap and power analysis tools
- Hierarchical/multilevel network support (HTNA/MLNA)
- 180+ learning action verbs for educational simulations

## Installation

```r
# Install from GitHub
devtools::install_github("mohsaqr/Saqrlab")
```

### Dependencies

```r
install.packages(c("tna", "seqHMM", "dplyr", "tidyr", "parallel",
                   "future", "future.apply", "progressr", "lhs"))
```

## Quick Start

```r
library(Saqrlab)
library(tna)

# Generate a TNA network in one line
model <- simulate_tna_network(seed = 42)

# Use with all tna functions
plot(model)
centralities(model)
communities(model)
```

## Function Reference

### Data Simulation

| Function | Description | Example |
|----------|-------------|---------|
| `simulate_tna_network()` | **Single fitted TNA model** | `simulate_tna_network(seed = 42)` |
| `simulate_matrix()` | Simple transition matrix | `simulate_matrix(n_nodes = 5, seed = 42)` |
| `simulate_htna()` | Multi-type HTNA/MLNA matrix | `simulate_htna(n_nodes = 5, n_types = 3)` |
| `simulate_sequences()` | Markov chain sequences | `simulate_sequences(n_sequences = 100)` |
| `simulate_sequences_advanced()` | Sequences with stability patterns | `simulate_sequences_advanced(...)` |
| `simulate_long_data()` | Hierarchical group data | `simulate_long_data(n_groups = 5)` |
| `simulate_onehot_data()` | One-hot encoded sequences | `simulate_onehot_data(n_sequences = 50)` |
| `simulate_edge_list()` | Social network edge lists | `simulate_edge_list(n_nodes = 20)` |

### Network Generation

| Function | Description | Example |
|----------|-------------|---------|
| `simulate_tna_datasets()` | Complete TNA datasets | `simulate_tna_datasets(n_datasets = 5)` |
| `simulate_tna_networks()` | Fitted TNA models | `simulate_tna_networks(n_networks = 5)` |
| `simulate_group_tna_networks()` | Group TNA models | `simulate_group_tna_networks(n_groups = 5)` |
| `simulate_tna_matrix()` | HTNA/MLNA matrices | `simulate_tna_matrix(n_states = 15)` |
| `generate_probabilities()` | Random transition probs | `generate_probabilities(n_states = 5)` |

### Model Fitting & Extraction

| Function | Description | Example |
|----------|-------------|---------|
| `fit_network_model()` | Fit TNA models | `fit_network_model(data, "tna")` |
| `extract_transition_matrix()` | Get transition matrix | `extract_transition_matrix(model)` |
| `extract_initial_probs()` | Get initial probabilities | `extract_initial_probs(model)` |
| `extract_edges()` | Get edge list | `extract_edges(model, threshold = 0.05)` |

### Network Comparison

| Function | Description | Example |
|----------|-------------|---------|
| `compare_networks()` | Compare two networks | `compare_networks(model1, model2)` |
| `compare_centralities()` | Compare centrality profiles | `compare_centralities(model1, model2)` |
| `compare_edge_recovery()` | Edge recovery metrics | `compare_edge_recovery(orig, sim)` |

### Data Conversion

| Function | Description | Example |
|----------|-------------|---------|
| `wide_to_long()` | Wide to long format | `wide_to_long(sequences)` |
| `long_to_wide()` | Long to wide format | `long_to_wide(data, "id", "Time", "Action")` |
| `prepare_for_tna()` | Prepare for tna package | `prepare_for_tna(data)` |
| `action_to_onehot()` | Convert to one-hot | `action_to_onehot(sequences)` |

### Batch Processing

| Function | Description | Example |
|----------|-------------|---------|
| `batch_fit_models()` | Fit multiple models | `batch_fit_models(data_list, "tna")` |
| `batch_apply()` | Apply function to list | `batch_apply(models, extract_edges)` |

### Bootstrap & Simulation

| Function | Description | Example |
|----------|-------------|---------|
| `run_bootstrap_simulation()` | Bootstrap analysis | `run_bootstrap_simulation(data, 100)` |
| `run_grid_simulation()` | Parameter grid search | `run_grid_simulation(param_grid)` |
| `run_network_simulation()` | Model comparison study | `run_network_simulation(...)` |
| `run_bootstrap_iteration()` | Evaluate bootstrap results | `run_bootstrap_iteration(results, model)` |
| `summarize_grid_results()` | Analyze grid output | `summarize_grid_results(grid_results)` |

### Learning States

| Function | Description | Example |
|----------|-------------|---------|
| `get_learning_states()` | Get verbs by category | `get_learning_states("metacognitive")` |
| `list_learning_categories()` | Show categories | `list_learning_categories()` |
| `select_states()` | Intelligent selection | `select_states(10)` |
| `LEARNING_STATES` | Full dataset (180+ verbs) | `LEARNING_STATES$cognitive` |
| `GLOBAL_NAMES` | 300 diverse names | `head(GLOBAL_NAMES, 20)` |

### Utilities

| Function | Description | Example |
|----------|-------------|---------|
| `generate_param_grid()` | Parameter grids | `generate_param_grid(n = c(50, 100))` |
| `validate_sim_params()` | Validate parameters | `validate_sim_params(params)` |
| `summarize_simulation()` | Summary stats | `summarize_simulation(results)` |
| `summarize_networks()` | Network summaries | `summarize_networks(model_list)` |

## Extensive Examples

### Example 1: Basic TNA Workflow

```r
library(Saqrlab)
library(tna)

# Simulate sequences with metacognitive states
sequences <- simulate_sequences(
  n_sequences = 200,
  seq_length = 25,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

# Fit TNA model
model <- fit_network_model(sequences, "tna")

# Analyze
plot(model)
centralities(model)

# Extract components
trans_mat <- extract_transition_matrix(model)
edges <- extract_edges(model, threshold = 0.05)
```

### Example 2: Hierarchical TNA (HTNA)

```r
library(Saqrlab)
library(tna)

# Generate 25-node multi-type matrix (5 types x 5 nodes)
net <- simulate_htna(seed = 42)

# View structure
dim(net$matrix)        # 25 x 25
names(net$node_types)  # 5 types

# Visualize with tna package
plot_htna(net$matrix, net$node_types, layout = "polygon")

# Custom configuration
net <- simulate_htna(
  n_nodes = 4,
  n_types = 3,
  type_names = c("Planning", "Execution", "Reflection"),
  within_prob = 0.5,    # High within-type connectivity
  between_prob = 0.15,  # Low between-type connectivity
  seed = 42
)
```

### Example 3: Group TNA Analysis

```r
library(Saqrlab)
library(tna)

# Generate hierarchical data (actors in groups in courses)
long_data <- simulate_long_data(
  n_groups = 5,
  n_actors = c(8, 12),    # Variable group sizes
  n_courses = 3,
  categories = "group_regulation",
  seq_length_range = c(10, 25),
  achiever_levels = c("High", "Low"),
  seed = 42
)

# View structure
head(long_data)
table(long_data$Group, long_data$Achiever)

# Convert to wide format
wide_data <- long_to_wide(
  long_data,
  id_col = "Actor",
  time_col = "Time",
  action_col = "Action"
)

# Fit model
model <- fit_network_model(wide_data, "tna")
```

### Example 4: Bootstrap Power Analysis

```r
library(Saqrlab)
library(tna)

# Create parameter grid for power analysis
power_grid <- generate_param_grid(
  param_ranges = list(
    n_sequences = c(25, 300),
    seq_length = c(15, 40),
    n_states = c(4, 8)
  ),
  n = 50,
  method = "lhs"
)

# Run power analysis (50 replications per condition)
power_results <- run_grid_simulation(
  param_grid = power_grid,
  n_runs_per_setting = 50,
  model_type = "tna",
  parallel = TRUE,
  seed = 42
)

# Analyze results
analysis <- summarize_grid_results(power_results)
summary_by_n <- summarize_simulation(
  power_results,
  by = "n_sequences",
  metrics = c("mean", "sd", "ci")
)

# Find minimum sample size for r >= 0.90
min_n <- summary_by_n$n_sequences[
  which(summary_by_n$correlation_mean >= 0.90)[1]
]
```

### Example 5: Model Comparison Study

```r
library(Saqrlab)
library(tna)

# Generate ground truth
true_probs <- generate_probabilities(n_states = 6, seed = 42)

# Simulate sequences
sequences <- simulate_sequences(
  trans_matrix = true_probs$transition_matrix,
  init_probs = true_probs$initial_probs,
  n_sequences = 200,
  seq_length = 30
)

# Fit all model types
models <- list(
  tna = fit_network_model(sequences, "tna"),
  ftna = fit_network_model(sequences, "ftna"),
  ctna = fit_network_model(sequences, "ctna"),
  atna = fit_network_model(sequences, "atna")
)

# Compare to TNA as reference
comparisons <- lapply(names(models)[-1], function(m) {
  comp <- compare_networks(models$tna, models[[m]])
  data.frame(
    model = m,
    correlation = comp$metrics$correlation,
    rmse = comp$metrics$rmse
  )
})
do.call(rbind, comparisons)

# Edge recovery analysis
recovery <- compare_edge_recovery(
  original = list(weights = true_probs$transition_matrix),
  simulated = models$tna,
  threshold = 0.05
)
cat("F1 Score:", recovery$f1_score)
```

### Example 6: Custom Learning States

```r
library(Saqrlab)

# View available categories
list_learning_categories()

# Get metacognitive verbs
meta_verbs <- get_learning_states("metacognitive")

# Get 10 random verbs from cognitive + behavioral
mixed_verbs <- get_learning_states(
  c("cognitive", "behavioral"),
  n = 10,
  seed = 42
)

# Smart selection based on network size
states_5 <- select_states(5, seed = 42)   # Single category
states_15 <- select_states(15, seed = 42) # Multiple categories

# Biased selection
srl_states <- select_states(
  n_states = 10,
  primary_categories = "metacognitive",
  secondary_categories = "cognitive",
  primary_ratio = 0.7,  # 70% metacognitive
  seed = 42
)

# Use in simulation
sequences <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 6,
  categories = c("metacognitive", "motivational"),
  seed = 42
)
```

## Learning States Reference

Saqrlab includes 180+ learning action verbs in 8 categories:

| Category | Count | Examples |
|----------|-------|----------|
| metacognitive | 20 | Plan, Monitor, Evaluate, Reflect, Regulate, Adjust |
| cognitive | 30 | Read, Study, Analyze, Summarize, Memorize, Connect |
| behavioral | 30 | Practice, Annotate, Research, Review, Revise, Write |
| social | 30 | Collaborate, Discuss, Seek_help, Question, Explain, Share |
| motivational | 30 | Focus, Persist, Explore, Create, Strive, Commit |
| affective | 30 | Enjoy, Appreciate, Value, Interest, Curious, Cope |
| group_regulation | 9 | Adapt, Cohesion, Consensus, Coregulate, Plan, Monitor |
| lms | 30 | View, Access, Download, Upload, Submit, Navigate |

```r
# Access directly
LEARNING_STATES$metacognitive

# Or use helper function
get_learning_states("social", n = 5, seed = 42)
```

## Documentation

- `?Saqrlab` - Package overview
- `vignette("introduction")` - Getting started
- `vignette("simulation-guide")` - Complete simulation guide
- `vignette("tna-workflow")` - End-to-end analysis
- `vignette("htna-mlna")` - Hierarchical/multilevel networks
- `vignette("bootstrap-power")` - Power analysis
- `vignette("learning-states")` - Learning states reference

## Citation

If you use Saqrlab in your research, please cite:

```
Saqr, M. (2025). Saqrlab: Simulation and Analysis Tools for Temporal Network
Analysis. R package version 0.1.0. https://github.com/mohsaqr/Saqrlab
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests at [GitHub](https://github.com/mohsaqr/Saqrlab).

## Author

**Mohammed Saqr** - [saqr@saqr.me](mailto:saqr@saqr.me)
