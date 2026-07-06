# Introduction to Saqrlab

## Overview

Saqrlab is an R package designed for simulating and analyzing Temporal
Network Analysis (TNA) models. It provides tools for:

- **Simulating** Markov chain sequences and transition matrices
- **Generating** complete TNA datasets with realistic parameters
- **Fitting** various TNA model types (TNA, fTNA, cTNA, aTNA)
- **Comparing** networks using multiple metrics
- **Running** bootstrap and simulation studies

The package is particularly useful for educational researchers working
with learning analytics data, as it includes built-in learning state
vocabularies for realistic simulations.

## Installation

``` r

# Install from GitHub
devtools::install_github("mohsaqr/Saqrlab")
```

## Getting Started

``` r

library(Saqrlab)
#> 
#> Attaching package: 'Saqrlab'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
```

### Your First Simulation

The simplest way to get started is to simulate sequences using the
built-in learning states:

``` r

# Simulate 100 sequences, each with 20 actions across 5 states
sequences <- simulate_sequences(
  n_sequences = 100,
  seq_length = 20,
  n_states = 5,
  seed = 42
)

# View the first few sequences
head(sequences)
#>           V1         V2         V3         V4         V5         V6         V7
#> 1 Appreciate      Judge Discourage Discourage   Regulate Discourage   Regulate
#> 2      Judge   Regulate      Judge       Plan Discourage Discourage Appreciate
#> 3      Judge      Judge   Regulate      Judge   Regulate Discourage Discourage
#> 4       Plan Discourage   Regulate Discourage      Judge   Regulate Appreciate
#> 5      Judge      Judge       Plan Appreciate   Regulate Discourage      Judge
#> 6 Discourage   Regulate      Judge Discourage      Judge   Regulate Appreciate
#>           V8         V9        V10        V11        V12        V13        V14
#> 1 Appreciate Discourage   Regulate Appreciate Discourage      Judge   Regulate
#> 2 Discourage   Regulate      Judge      Judge Discourage Discourage      Judge
#> 3      Judge      Judge   Regulate Appreciate Discourage      Judge   Regulate
#> 4   Regulate Discourage   Regulate Discourage      Judge      Judge Discourage
#> 5      Judge   Regulate Discourage      Judge   Regulate Appreciate   Regulate
#> 6 Discourage      Judge   Regulate Discourage Discourage   Regulate Discourage
#>          V15        V16        V17        V18        V19        V20
#> 1 Appreciate Discourage   Regulate Discourage   Regulate Appreciate
#> 2      Judge Discourage      Judge   Regulate Discourage Discourage
#> 3 Discourage Discourage Appreciate Discourage   Regulate Discourage
#> 4 Discourage Discourage Appreciate      Judge   Regulate Discourage
#> 5 Appreciate Discourage Discourage   Regulate Appreciate Discourage
#> 6   Regulate Appreciate Discourage      Judge Discourage Appreciate

# See the state names (randomly selected learning verbs)
colnames(unique(as.matrix(sequences)))
#>  [1] "V1"  "V2"  "V3"  "V4"  "V5"  "V6"  "V7"  "V8"  "V9"  "V10" "V11" "V12"
#> [13] "V13" "V14" "V15" "V16" "V17" "V18" "V19" "V20"
```

### Fitting a TNA Model

Once you have sequences, you can fit a TNA model:

``` r

library(tna)

# Fit a TNA model
model <- fit_network_model(sequences, "tna")

# View the model
print(model)

# Plot the network
plot(model)
```

### Extracting Network Components

You can extract various components from fitted models:

``` r

# Get the transition matrix
trans_mat <- extract_transition_matrix(model)
print(trans_mat)

# Get initial probabilities
init_probs <- extract_initial_probs(model)
print(init_probs)

# Get edge list
edges <- extract_edges(model, threshold = 0.05)
head(edges)
```

## Key Concepts

### Learning States

Saqrlab includes over 180 learning action verbs organized into 8
categories:

``` r

# View available categories
list_learning_categories()
#>           Category Count
#> 1    metacognitive    20
#> 2        cognitive    30
#> 3       behavioral    30
#> 4           social    30
#> 5     motivational    30
#> 6        affective    30
#> 7 group_regulation     9
#> 8              lms    30
#>                                                  Examples
#> 1         Plan, Monitor, Evaluate, Reflect, Regulate, ...
#> 2          Read, Study, Analyze, Summarize, Memorize, ...
#> 3       Practice, Annotate, Research, Review, Revise, ...
#> 4 Collaborate, Discuss, Seek_help, Question, Explain, ...
#> 5            Focus, Persist, Explore, Create, Strive, ...
#> 6        Enjoy, Appreciate, Value, Interest, Curious, ...
#> 7    Adapt, Cohesion, Consensus, Coregulate, Discuss, ...
#> 8             View, Access, Download, Upload, Submit, ...
```

You can retrieve verbs from specific categories:

``` r

# Get metacognitive verbs
get_learning_states("metacognitive")
#>  [1] "Plan"        "Monitor"     "Evaluate"    "Reflect"     "Regulate"   
#>  [6] "Adjust"      "Adapt"       "Check"       "Assess"      "Judge"      
#> [11] "Strategize"  "Prioritize"  "Set_goals"   "Track"       "Self_assess"
#> [16] "Calibrate"   "Diagnose"    "Forecast"    "Anticipate"  "Reconsider"

# Get 6 random verbs from cognitive and behavioral
get_learning_states(c("cognitive", "behavioral"), n = 6, seed = 42)
#> [1] "Submit"     "Write"      "Read"       "Generalize" "Compare"   
#> [6] "Test"
```

### Transition Matrices

A transition matrix defines the probability of moving from one state to
another. Each row sums to 1:

``` r

# Generate a simple transition matrix
mat <- simulate_matrix(n_nodes = 4, seed = 42)
print(mat)
#>          Regulate Plan Judge Reflect
#> Regulate        0    0     1       0
#> Plan            1    0     0       0
#> Judge           1    0     0       0
#> Reflect         1    0     0       0

# Verify rows sum to 1
rowSums(mat)
#> Regulate     Plan    Judge  Reflect 
#>        1        1        1        1
```

### Multi-Type Networks (HTNA/MLNA)

For hierarchical or multilevel network analysis, you can generate
matrices with multiple node types:

``` r

# Generate a 15-node matrix with 3 types (5 nodes each)
net <- simulate_htna(n_nodes = 5, n_types = 3, seed = 42)

# View the node types
net$node_types
#> $Metacognitive
#> [1] "Diagnose" "Regulate" "Plan"     "Judge"    "Reflect" 
#> 
#> $Cognitive
#> [1] "Understand" "Classify"   "Process"    "Encode"     "Abstract"  
#> 
#> $Behavioral
#> [1] "Write"   "Review"  "Outline" "Revise"  "Draft"
```

## Basic Workflows

### Workflow 1: Simulate and Analyze

``` r

library(Saqrlab)
library(tna)

# 1. Simulate sequences
sequences <- simulate_sequences(
  n_sequences = 200,
  seq_length = 25,
  n_states = 6,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)

# 2. Fit TNA model
model <- fit_network_model(sequences, "tna")

# 3. Analyze
plot(model)
centralities(model)
```

### Workflow 2: Compare Two Networks

``` r

# Create two different networks
net1 <- simulate_tna_networks(1, n_states = 5, seed = 42)
net2 <- simulate_tna_networks(1, n_states = 5, seed = 123)

# Compare them
comparison <- compare_networks(
  net1$network_1$model,
  net2$network_1$model
)

# View metrics
comparison$metrics
```

### Workflow 3: Bootstrap Analysis

``` r

# Run bootstrap simulation
results <- run_bootstrap_simulation(
  original_data = sequences,
  n_bootstrap = 50,
  model_type = "tna",
  seed = 42
)
```

## Data Formats

Saqrlab works with multiple data formats:

### Wide Format (Default)

Each row is a sequence, columns are time points:

       V1       V2        V3      V4
    1  Plan     Monitor   Read    Plan
    2  Read     Plan      Study   Monitor

### Long Format

One observation per row with ID, time, and action columns:

      id  Time  Action
      1   1     Plan
      1   2     Monitor
      1   3     Read

Convert between formats:

``` r

# Wide to long
long_data <- wide_to_long(sequences)

# Long to wide
wide_data <- long_to_wide(long_data, id_col = "id",
                          time_col = "Time", action_col = "Action")
```

## Next Steps

- **Simulation Guide**: Learn all simulation functions in detail
- **TNA Workflow**: Complete end-to-end analysis examples
- **HTNA & MLNA**: Hierarchical and multilevel network analysis
- **Bootstrap Power**: Power analysis for sample size planning
- **Learning States**: Full reference for educational research

## Getting Help

- Package documentation:
  [`?Saqrlab`](https://pak.dynasite.org/Saqrlab/reference/Saqrlab-package.md)
- Function help: `?function_name`
- GitHub issues: <https://github.com/mohsaqr/Saqrlab/issues>
