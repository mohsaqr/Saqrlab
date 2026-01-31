# Saqrlab

**Simulation and Analysis Tools for Temporal Network Analysis (TNA)**

[![R CMD check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Saqrlab is an R package that provides comprehensive tools for simulating Markov chain sequences, bootstrapping temporal network models, and evaluating network recovery performance. It is designed to work seamlessly with the [`tna`](https://github.com/sonsoleslp/tna) package for Temporal Network Analysis.

## Table of Contents

- [Installation](#installation)
- [Overview](#overview)
- [Quick Start](#quick-start)
- [Functions Reference](#functions-reference)
  - [Sequence Generation](#sequence-generation)
  - [Network Generation](#network-generation)
  - [Model Fitting](#model-fitting)
  - [Data Conversion](#data-conversion)
  - [Network Comparison](#network-comparison)
  - [Model Extraction](#model-extraction)
  - [Simulation & Bootstrap](#simulation--bootstrap)
  - [Batch Processing](#batch-processing)
  - [Summary Functions](#summary-functions)
  - [Learning States](#learning-states)
  - [Parameter Handling](#parameter-handling)
- [Detailed Examples](#detailed-examples)
- [Workflows](#workflows)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [License](#license)

## Installation

### From GitHub

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install Saqrlab from GitHub
devtools::install_github("mohsaqr/Saqrlab")
```
### Dependencies

Saqrlab requires the following packages:

```r
install.packages(c("tna", "seqHMM", "dplyr", "tidyr", "parallel",
                   "future", "future.apply", "progressr", "lhs"))
```

## Overview

Saqrlab provides a complete toolkit for TNA simulation studies:

| Category | Functions | Purpose |
|----------|-----------|---------|
| **Sequence Generation** | `simulate_sequences`, `simulate_sequences_advanced`, `simulate_long_data` | Generate Markov chain sequences in various formats |
| **Network Generation** | `generate_tna_networks`, `generate_probabilities` | Create TNA networks with realistic parameters |
| **Data Conversion** | `wide_to_long`, `long_to_wide`, `prepare_for_tna` | Convert between data formats |
| **Network Comparison** | `compare_networks`, `compare_centralities`, `calculate_edge_recovery` | Evaluate network similarity and recovery |
| **Model Extraction** | `extract_transition_matrix`, `extract_initial_probs`, `extract_edges` | Extract components from TNA models |
| **Simulation** | `run_network_simulation`, `run_bootstrap_simulation`, `run_grid_simulation` | Run comprehensive simulation studies |
| **Batch Processing** | `batch_fit_models`, `batch_apply` | Process multiple datasets efficiently |
| **Summary** | `summarize_simulation`, `summarize_networks` | Aggregate and summarize results |
| **Learning States** | `get_learning_states`, `smart_select_states`, `LEARNING_STATES` | Realistic educational state names |

## Quick Start

```r
library(Saqrlab)
library(tna)

# 1. Create a transition matrix
trans_mat <- matrix(c(
  0.7, 0.2, 0.1,
  0.3, 0.5, 0.2,
  0.2, 0.3, 0.5
), nrow = 3, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("A", "B", "C")

# 2. Define initial probabilities
init_probs <- c(A = 0.5, B = 0.3, C = 0.2)

# 3. Simulate sequences
sequences <- simulate_sequences(
 transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 20,
  num_rows = 100
)

# 4. Fit a TNA model
model <- fit_network_model(sequences, "tna")

# 5. View the network
plot(model)
```

## Functions Reference

### Sequence Generation

#### `simulate_sequences()`

Generate Markov chain sequences. Can use provided parameters or auto-generate random ones with learning state names.

```r
# Method 1: Provide your own transition matrix
sequences <- simulate_sequences(
  transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 20,
  num_rows = 100,
  min_na = 0,
  max_na = 5,        # Add 0-5 trailing NAs per sequence
  include_na = TRUE
)

# Method 2: Auto-generate with letter names (A, B, C, ...)
sequences <- simulate_sequences(
  max_seq_length = 20,
  num_rows = 100,
  n_states = 5,
  seed = 42
)

# Method 3: Auto-generate with learning state names
sequences <- simulate_sequences(
  max_seq_length = 25,
  num_rows = 150,
  n_states = 6,
  use_learning_states = TRUE,
  learning_categories = c("metacognitive", "cognitive"),
  seed = 123
)

head(sequences)
#>        V1      V2       V3     V4        V5 ...
#> 1    Plan Monitor     Read  Study  Evaluate ...
#> 2    Read    Plan  Monitor   Plan      Read ...
#> 3 Monitor    Read Evaluate   Read      Plan ...

# Method 4: Get sequences AND generating parameters back
result <- simulate_sequences(
  max_seq_length = 20,
  num_rows = 100,
  n_states = 4,
  use_learning_states = TRUE,
  return_params = TRUE,
  seed = 42
)
result$sequences          # The sequences
result$transition_matrix  # The random transition matrix used
result$state_names        # e.g., c("Plan", "Monitor", "Read", "Practice")

# Method 5: Custom state names
sequences <- simulate_sequences(
  max_seq_length = 20,
  num_rows = 100,
  n_states = 4,
  state_names = c("Explore", "Learn", "Practice", "Master"),
  seed = 42
)
```

#### `simulate_sequences_advanced()`

Generate sequences with stability modes for more realistic patterns.

```r
# Define stable transition patterns
stable_transitions <- list(
  c("Plan", "Monitor"),
  c("Monitor", "Evaluate"),
  c("Read", "Summarize")
)

# Simulate with 90% stability
sequences_adv <- simulate_sequences_advanced(
  transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 30,
  num_rows = 200,
  stable_transitions = stable_transitions,
  stability_prob = 0.90,
  unstable_mode = "random_jump",  # or "perturb_prob", "unlikely_jump"
  unstable_random_transition_prob = 0.3
)
```

#### `simulate_long_data()`

Generate long-format data with hierarchical structure (actors in groups in courses).

```r
# Simulate educational data with group structure
long_data <- simulate_long_data(
  n_groups = 20,
  actors_per_group = 10,       # or c(8, 12) for variable sizes
  n_courses = 3,
  actions = NULL,              # Use default learning actions
  action_categories = "group_regulation",
  seq_length_range = c(8, 25),
  achiever_levels = c("High", "Low"),
  achiever_probs = c(0.5, 0.5),
  seed = 42
)

head(long_data)
#>   Actor Achiever Group Course                Time   Action
#> 1     1     High     1      A 2025-01-01 10:05:23    plan
#> 2     1     High     1      A 2025-01-01 10:12:45  monitor
#> 3     1     High     1      A 2025-01-01 10:18:12  discuss
```

### Network Generation

#### `generate_tna_networks()`

Generate multiple TNA networks with random or specified parameters.

```r
# Generate 5 networks with letter state names
networks <- generate_tna_networks(
  n_networks = 5,
  n_states = 4,
  num_rows = 100,
  max_seq_length = 30,
  model_type = "tna",
  seed = 42,
  verbose = TRUE
)

# Access network components
networks$network_1$model              # The TNA model
networks$network_1$transition_probs   # Generating transition matrix
networks$network_1$initial_probs      # Generating initial probabilities
networks$network_1$sequences          # The simulated sequences

# Generate networks with realistic learning state names
learning_networks <- generate_tna_networks(
  n_networks = 5,
  n_states = 8,
  use_learning_states = TRUE,
  learning_categories = c("metacognitive", "cognitive"),
  smart_select = TRUE,
  seed = 123
)

# View the learning state names used
learning_networks$network_1$params$state_names
#> [1] "Plan" "Monitor" "Evaluate" "Read" "Analyze" "Summarize" "Reflect" "Study"
```

#### `generate_probabilities()`

Generate random transition matrices and initial probabilities.
```r
probs <- generate_probabilities(n_states = 5)
probs$transition_matrix
probs$initial_probs
```

### Model Fitting

#### `fit_network_model()`

Fit TNA models to sequence data.

```r
# Fit different model types
model_tna  <- fit_network_model(sequences, "tna")   # Standard TNA
model_ftna <- fit_network_model(sequences, "ftna")  # Filtered TNA
model_ctna <- fit_network_model(sequences, "ctna")  # Conditional TNA
model_atna <- fit_network_model(sequences, "atna")  # Aggregated TNA
```

### Data Conversion

#### `wide_to_long()`

Convert wide-format sequences to long format.

```r
# Convert sequences to long format
long_data <- wide_to_long(
  data = sequences,
  id_col = NULL,          # Auto-generate IDs
  time_prefix = "V",      # Column prefix (V1, V2, ...)
  action_col = "Action",
  time_col = "Time",
  drop_na = TRUE
)

head(long_data)
#>   id Time Action
#> 1  1    1      A
#> 2  1    2      A
#> 3  1    3      B
```

#### `long_to_wide()`

Convert long-format data to wide format.

```r
# Convert back to wide format
wide_data <- long_to_wide(
  data = long_data,
  id_col = "id",
  time_col = "Time",
  action_col = "Action",
  time_prefix = "V",
  fill_na = TRUE
)
```

#### `prepare_for_tna()`

Prepare any data format for TNA analysis.

```r
# Auto-detect format and prepare for TNA
tna_ready <- prepare_for_tna(
  data = my_data,
  type = "auto",          # or "sequences", "long"
  state_names = c("A", "B", "C"),
  validate = TRUE
)

model <- tna::tna(tna_ready)
```

### Network Comparison

#### `compare_networks()`

Compare two TNA networks using multiple metrics.

```r
# Compare original and simulated models
comparison <- compare_networks(
  model1 = original_model,
  model2 = simulated_model,
  metrics = c("correlation", "rmse", "mae", "edge_diff", "cosine"),
  scaling = "none",       # or "minmax", "zscore"
  include_self = TRUE,
  threshold = 0.05
)

# View metrics
comparison$metrics
#> $correlation
#> [1] 0.9234
#> $rmse
#> [1] 0.0456
#> $edge_diff
#> [1] 0.1111

# View edge-level comparison
head(comparison$edge_comparison)
#>   from to    weight1    weight2         diff     abs_diff
#> 1    A  A 0.70000000 0.68234521 -0.017654789 0.0176547890
#> 2    B  A 0.30000000 0.31234567  0.012345670 0.0123456700

# Print summary
cat(comparison$summary)
```

#### `compare_centralities()`

Compare centrality profiles between networks.

```r
cent_comparison <- compare_centralities(
  model1 = original_model,
  model2 = simulated_model,
  measures = c("OutStrength", "InStrength", "Betweenness"),
  method = "both"  # Pearson and Spearman correlations
)

cent_comparison$correlations
#> $OutStrength_pearson
#> [1] 0.9512
#> $OutStrength_spearman
#> [1] 0.9000
#> $InStrength_pearson
#> [1] 0.8934
```

#### `calculate_edge_recovery()`

Calculate edge recovery metrics (precision, recall, F1).

```r
recovery <- calculate_edge_recovery(
  original = original_model,
  simulated = simulated_model,
  threshold = 0.01,
  return_edges = TRUE
)

# Classification metrics
recovery$precision   # TP / (TP + FP)
recovery$recall      # TP / (TP + FN)
recovery$f1_score    # Harmonic mean
recovery$accuracy    # (TP + TN) / total
recovery$jaccard     # TP / (TP + FP + FN)

# Confusion matrix counts
recovery$true_positives
recovery$false_positives
recovery$false_negatives
recovery$true_negatives

# Edge-level details
head(recovery$edges)
#>   from to weight_original weight_simulated present_original present_simulated status
#> 1    A  A      0.70000000       0.68234521             TRUE              TRUE     TP
#> 2    B  A      0.30000000       0.31234567             TRUE              TRUE     TP
#> 3    C  A      0.00500000       0.00000000            FALSE             FALSE     TN
```

### Model Extraction

#### `extract_transition_matrix()`

Extract transition matrix from a TNA model.

```r
# Get raw weights
trans_mat <- extract_transition_matrix(model, type = "raw")

# Get row-normalized (rows sum to 1)
trans_mat_scaled <- extract_transition_matrix(model, type = "scaled")
rowSums(trans_mat_scaled)  # All 1s
```

#### `extract_initial_probs()`

Extract initial state probabilities.

```r
init_probs <- extract_initial_probs(model)
#>         A         B         C
#> 0.4523456 0.3234567 0.2241977
```

#### `extract_edges()`

Extract edge list from a model.

```r
edges <- extract_edges(
  model,
  threshold = 0.05,      # Minimum weight
  include_self = FALSE,  # Exclude self-loops
  sort_by = "weight"     # Sort by weight (descending)
)

head(edges)
#>   from to    weight
#> 1    A  A 0.7000000
#> 2    B  B 0.5000000
#> 3    C  C 0.5000000
#> 4    A  B 0.2000000
```

### Simulation & Bootstrap

#### `run_network_simulation()`

Run comprehensive network simulations.

```r
# Run simulations comparing model types
results <- run_network_simulation(
  original_data_list = list(data1, data2, data3),
  sim_params = list(
    list(max_seq_length = 20, num_rows = 50),
    list(max_seq_length = 30, num_rows = 100),
    list(max_seq_length = 40, num_rows = 200)
  ),
  models = c("tna", "ftna"),
  comparisons = c("original", "across_models"),
  num_runs = 10,
  parallel = TRUE,
  scaling = "minmax"
)

# View results
results$metrics       # All comparison metrics
results$summary_stats # Aggregated statistics
```

#### `run_bootstrap_simulation()`

Run bootstrap simulations for stability analysis.

```r
bootstrap_results <- run_bootstrap_simulation(
  original_data = sequences,
  n_bootstrap = 100,
  model_type = "tna",
  sample_size = NULL,  # Use original size
  parallel = TRUE,
  seed = 42
)
```

#### `run_grid_simulation()`

Run simulations across a parameter grid.

```r
# Create parameter grid
param_grid <- create_param_grid(
  num_rows = c(50, 100, 200),
  max_seq_length = c(15, 30),
  num_states = c(4, 6, 8)
)

# Run grid simulation
grid_results <- run_grid_simulation(
  param_grid = param_grid,
  n_runs_per_setting = 20,
  model_type = "tna",
  parallel = TRUE
)

# Analyze results
analysis <- analyze_grid_results(grid_results)
```

### Batch Processing

#### `batch_fit_models()`

Fit models to multiple datasets.

```r
# Generate multiple datasets
datasets <- lapply(1:20, function(i) {
  simulate_sequences(trans_mat, init_probs, 20, 100)
})

# Fit models in parallel
models <- batch_fit_models(
  data_list = datasets,
  model_type = "tna",
  parallel = TRUE,
  cores = 4,
  progress = TRUE
)
#> Fitting 20 tna models using 4 cores...
#> Successfully fitted 20/20 models
```

#### `batch_apply()`

Apply any function to multiple models.

```r
# Extract all transition matrices
trans_matrices <- batch_apply(
  object_list = models,
  fun = extract_transition_matrix,
  parallel = TRUE
)

# Compare each model to a reference
ref_model <- models[[1]]
correlations <- batch_apply(
  models[-1],
  function(m) compare_networks(ref_model, m)$metrics$correlation,
  simplify = TRUE
)
mean(correlations)
```

### Summary Functions

#### `summarize_simulation()`

Summarize simulation results.

```r
# Overall summary
summary_all <- summarize_simulation(
  results = simulation_results,
  metrics = c("mean", "sd", "ci")
)

# Summary by parameter
summary_by_n <- summarize_simulation(
  results = simulation_results,
  by = "num_rows",
  metrics = "all"  # mean, sd, median, ci, min, max, n
)

summary_by_n
#>   num_rows correlation_mean correlation_sd correlation_ci_lower ...
#> 1       50           0.8234         0.0456              0.8023
#> 2      100           0.9123         0.0234              0.9012
#> 3      200           0.9567         0.0123              0.9489
```

#### `summarize_networks()`

Summarize multiple network models.

```r
network_summary <- summarize_networks(
  model_list = models,
  include = c("density", "centrality", "edges"),
  threshold = 0.01,
  centrality_measures = c("OutStrength", "InStrength")
)

# Per-network statistics
network_summary$summary_table
#>   network_id   density edge_mean   edge_sd edge_max OutStrength_mean ...
#> 1          1 0.8888889 0.1111111 0.1847521 0.680234        0.3333333
#> 2          2 0.9444444 0.1111111 0.1623456 0.712345        0.3333333

# Aggregate statistics
network_summary$aggregate$density
#> $mean [1] 0.9166667
#> $sd   [1] 0.0392837
#> $median [1] 0.9166667
```

### Learning States

Saqrlab includes 180+ realistic learning action verbs organized into 7 categories for educational research simulations.

#### Categories

| Category | Count | Examples |
|----------|-------|----------|
| metacognitive | 20 | Plan, Monitor, Evaluate, Reflect, Regulate |
| cognitive | 30 | Read, Study, Analyze, Summarize, Memorize |
| behavioral | 30 | Practice, Annotate, Research, Review, Write |
| social | 30 | Collaborate, Discuss, Explain, Share, Teach |
| motivational | 30 | Focus, Persist, Explore, Strive, Commit |
| affective | 30 | Enjoy, Appreciate, Cope, Curious, Manage |
| group_regulation | 9 | Adapt, Cohesion, Consensus, Coregulate, Plan |

#### `get_learning_states()`

```r
# Get all metacognitive verbs
get_learning_states("metacognitive")

# Get verbs from multiple categories
get_learning_states(c("cognitive", "behavioral"))

# Random sample of 10 verbs
get_learning_states(n = 10, seed = 42)

# 8 verbs from social and motivational
get_learning_states(c("social", "motivational"), n = 8)
```

#### `smart_select_states()`

Intelligently select states based on network size.

```r
# 5 states (auto-selects 1 category)
smart_select_states(5, seed = 42)

# 10 states focused on metacognitive
smart_select_states(10,
  primary_categories = "metacognitive",
  secondary_categories = "cognitive",
  primary_ratio = 0.6
)

# 20 states balanced across all categories
smart_select_states(20, seed = 123)
```

#### `list_learning_categories()`

```r
list_learning_categories()
#>          Category Count                           Examples
#> 1   metacognitive    20 Plan, Monitor, Evaluate, Reflect, Regulate, ...
#> 2       cognitive    30 Read, Study, Analyze, Summarize, Memorize, ...
#> 3      behavioral    30 Practice, Annotate, Research, Review, Revise, ...
#> 4          social    30 Collaborate, Discuss, Seek_help, Question, Explain, ...
#> 5    motivational    30 Focus, Persist, Explore, Create, Strive, ...
#> 6       affective    30 Enjoy, Appreciate, Value, Interest, Curious, ...
#> 7 group_regulation     9 Adapt, Cohesion, Consensus, Coregulate, Discuss, ...
```

### Parameter Handling

#### `validate_sim_params()`

Validate and set defaults for simulation parameters.

```r
params <- validate_sim_params(list(
  num_rows = 100,
  max_seq_length = 30
))
# Adds defaults for min_na, max_na, etc.
```

#### `create_param_grid()`

Create parameter combinations for grid simulations.

```r
grid <- create_param_grid(
  num_rows = c(50, 100, 200),
  max_seq_length = c(15, 30, 45),
  num_states = c(4, 6)
)

nrow(grid)  # 18 combinations
```

## Detailed Examples

### Example 1: Complete Simulation Study

```r
library(Saqrlab)
library(tna)

# 1. Generate ground truth networks with learning states
set.seed(42)
networks <- generate_tna_networks(
  n_networks = 10,
  n_states = 6,
  use_learning_states = TRUE,
  learning_categories = c("metacognitive", "cognitive"),
  num_rows = 200,
  max_seq_length = 30,
  model_type = "tna",
  include_probabilities = TRUE
)

# 2. For each network, simulate recovery at different sample sizes
results <- lapply(networks, function(net) {
  original_trans <- net$transition_probs
  original_init <- net$initial_probs

  # Test different sample sizes
  sample_sizes <- c(50, 100, 200, 500)

  recovery_results <- lapply(sample_sizes, function(n) {
    # Simulate new sequences
    new_seqs <- simulate_sequences(
      transition_matrix = original_trans,
      initial_probabilities = original_init,
      max_seq_length = 30,
      num_rows = n
    )

    # Fit model
    new_model <- fit_network_model(new_seqs, "tna")

    # Compare to original
    comparison <- compare_networks(net$model, new_model)
    recovery <- calculate_edge_recovery(net$model, new_model)

    data.frame(
      sample_size = n,
      correlation = comparison$metrics$correlation,
      rmse = comparison$metrics$rmse,
      precision = recovery$precision,
      recall = recovery$recall,
      f1 = recovery$f1_score
    )
  })

  do.call(rbind, recovery_results)
})

# 3. Combine and summarize
all_results <- do.call(rbind, results)
summary_by_size <- summarize_simulation(
  all_results,
  by = "sample_size",
  metrics = c("mean", "sd", "ci")
)
print(summary_by_size)
```

### Example 2: Bootstrap Analysis

```r
library(Saqrlab)
library(tna)

# Generate original data
trans_mat <- matrix(c(
  0.6, 0.3, 0.1,
  0.2, 0.5, 0.3,
  0.1, 0.2, 0.7
), nrow = 3, byrow = TRUE)
rownames(trans_mat) <- colnames(trans_mat) <- c("Plan", "Execute", "Review")
init_probs <- c(Plan = 0.5, Execute = 0.3, Review = 0.2)

original_data <- simulate_sequences(
  transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 25,
  num_rows = 150
)

# Fit original model
original_model <- fit_network_model(original_data, "tna")

# Run bootstrap
bootstrap_results <- run_bootstrap_simulation(
  original_data = original_data,
  n_bootstrap = 100,
  model_type = "tna",
  parallel = TRUE,
  seed = 123
)

# Evaluate bootstrap performance
evaluation <- evaluate_bootstrap(
  bootstrap_results = bootstrap_results,
  original_model = original_model
)

print(evaluation)
```

### Example 3: Comparing Model Types

```r
library(Saqrlab)
library(tna)

# Generate sequences
sequences <- simulate_sequences(
  transition_matrix = trans_mat,
  initial_probabilities = init_probs,
  max_seq_length = 30,
  num_rows = 200
)

# Fit all model types
models <- list(
  tna = fit_network_model(sequences, "tna"),
  ftna = fit_network_model(sequences, "ftna"),
  ctna = fit_network_model(sequences, "ctna"),
  atna = fit_network_model(sequences, "atna")
)

# Compare all models to TNA as reference
comparisons <- lapply(models[-1], function(m) {
  compare_networks(models$tna, m)
})

# Extract correlations
sapply(comparisons, function(c) c$metrics$correlation)
#>      ftna      ctna      atna
#> 0.9876543 0.8765432 0.9234567
```

### Example 4: Working with Long Format Data

```r
library(Saqrlab)
library(tna)

# Simulate educational data
edu_data <- simulate_long_data(
  n_groups = 50,
  actors_per_group = c(8, 12),  # Variable group sizes
  n_courses = 3,
  action_categories = "group_regulation",
  seq_length_range = c(10, 30),
  seed = 42
)

# Check structure
str(edu_data)
table(edu_data$Course)
table(edu_data$Achiever)

# Convert to wide format for TNA
wide_data <- long_to_wide(
  data = edu_data,
  id_col = "Actor",
  time_col = "Time",
  action_col = "Action"
)

# Fit model
model <- fit_network_model(wide_data, "tna")
plot(model)

# Or prepare automatically
tna_data <- prepare_for_tna(edu_data, type = "long")
model2 <- tna::tna(tna_data)
```

## Workflows

### Workflow 1: Power Analysis for TNA Studies

Determine how many sequences you need for reliable network recovery.

```r
# Create parameter grid for power analysis
power_grid <- create_param_grid(
  num_rows = seq(25, 300, by = 25),
  max_seq_length = c(15, 25, 40),
  num_states = c(4, 6, 8)
)

# Run simulations
power_results <- run_grid_simulation(
  param_grid = power_grid,
  n_runs_per_setting = 50,
  model_type = "tna",
  parallel = TRUE
)

# Analyze to find minimum sample size for 0.90 correlation
analysis <- analyze_grid_results(power_results)
```

### Workflow 2: Model Comparison Study

Compare performance of different TNA model types.

```r
# Generate diverse test networks
test_networks <- generate_tna_networks(
  n_networks = 20,
  n_states = 6,
  use_learning_states = TRUE,
  num_rows = 150,
  seed = 42
)

# Fit all model types to each network's sequences
model_comparison <- lapply(test_networks, function(net) {
  seqs <- net$sequences

  models <- list(
    tna = fit_network_model(seqs, "tna"),
    ftna = fit_network_model(seqs, "ftna")
  )

  # Compare each to generating parameters
  original <- list(weights = net$transition_probs)

  list(
    tna_recovery = calculate_edge_recovery(original, models$tna),
    ftna_recovery = calculate_edge_recovery(original, models$ftna)
  )
})
```

## Dependencies

- **tna**: Temporal Network Analysis models
- **seqHMM**: Sequence analysis and HMM utilities
- **dplyr**: Data manipulation
- **tidyr**: Data tidying
- **parallel**: Parallel processing
- **future**: Asynchronous evaluation
- **future.apply**: Apply functions with futures
- **progressr**: Progress reporting
- **lhs**: Latin hypercube sampling for parameter grids

## Citation

If you use Saqrlab in your research, please cite:

```
Saqr, M. (2025). Saqrlab: Simulation and Analysis Tools for Temporal Network
Analysis. R package version 0.1.0. https://github.com/mohsaqr/Saqrlab
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Author

**Mohammed Saqr** - [saqr@saqr.me](mailto:saqr@saqr.me)
