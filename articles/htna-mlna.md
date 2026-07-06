# Hierarchical and Multilevel Network Analysis

``` r

library(Saqrlab)
#> 
#> Attaching package: 'Saqrlab'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
```

## Overview

This vignette covers hierarchical (HTNA) and multilevel (MLNA) network
analysis using Saqrlab. These approaches extend standard TNA by
organizing nodes into types or levels:

- **HTNA (Hierarchical TNA)**: Nodes grouped into distinct
  types/categories
- **MLNA (Multilevel TNA)**: Nodes organized into hierarchical levels
- **MTNA (Multi-Type TNA)**: General term for networks with typed nodes

Saqrlab provides
[`simulate_htna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md),
[`simulate_mlna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md),
and
[`simulate_mtna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
as aliases that produce identical output.

## Understanding Multi-Type Networks

### What Makes HTNA/MLNA Different?

In standard TNA, all nodes are treated equally. In HTNA/MLNA:

1.  **Nodes have types**: Each node belongs to a category
2.  **Within-type patterns**: Transitions within the same type
3.  **Between-type patterns**: Transitions across different types
4.  **Structured visualization**: Layout reflects node groupings

### Educational Research Applications

Multi-type networks are useful for studying:

- **Self-regulated learning**: Metacognitive, cognitive, behavioral
  dimensions
- **Collaborative learning**: Individual vs. group regulation
- **Multimodal data**: Different data sources (logs, surveys,
  observations)
- **Hierarchical processes**: Macro, meso, micro levels of analysis

## Simulating Multi-Type Networks

### Basic Usage

``` r

# Default: 5 types x 5 nodes = 25-node network
net <- simulate_htna(seed = 42)

# Network dimensions
dim(net$matrix)
#> [1] 25 25

# Node types
names(net$node_types)
#> [1] "Metacognitive" "Cognitive"     "Behavioral"    "Social"       
#> [5] "Motivational"

# Nodes per type
net$n_nodes_per_type
#> Metacognitive     Cognitive    Behavioral        Social  Motivational 
#>             5             5             5             5             5
```

### Understanding the Output

The function returns a list with:

``` r

# Transition matrix
head(net$matrix[1:6, 1:6])
#>            Diagnose Regulate Plan  Judge Reflect Understand
#> Diagnose     0.0000   0.0000    0 0.0000  0.0000     0.0000
#> Regulate     0.3553   0.0000    0 0.1239  0.0000     0.0000
#> Plan         0.0000   0.0000    0 0.0000  0.0742     0.2517
#> Judge        0.0000   0.0000    0 0.0000  0.0000     0.0000
#> Reflect      0.0000   0.4497    0 0.0472  0.0000     0.0000
#> Understand   0.0000   0.0000    0 0.0000  0.1433     0.0000

# Node types (for visualization)
net$node_types
#> $Metacognitive
#> [1] "Diagnose" "Regulate" "Plan"     "Judge"    "Reflect" 
#> 
#> $Cognitive
#> [1] "Understand" "Classify"   "Process"    "Encode"     "Abstract"  
#> 
#> $Behavioral
#> [1] "Write"   "Review"  "Outline" "Revise"  "Draft"  
#> 
#> $Social
#> [1] "Help"       "Critique"   "Contribute" "Present"    "Seek_help" 
#> 
#> $Motivational
#> [1] "Overcome"   "Accomplish" "Improve"    "Create"     "Strive"

# Type names
net$type_names
#> [1] "Metacognitive" "Cognitive"     "Behavioral"    "Social"       
#> [5] "Motivational"
```

### Default Learning Categories

By default,
[`simulate_htna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
uses five learning state categories:

| Type          | Learning States                         |
|---------------|-----------------------------------------|
| Metacognitive | Plan, Monitor, Evaluate, Reflect, …     |
| Cognitive     | Read, Study, Analyze, Summarize, …      |
| Behavioral    | Practice, Annotate, Research, Review, … |
| Social        | Collaborate, Discuss, Explain, Share, … |
| Motivational  | Focus, Persist, Explore, Create, …      |

## Customizing Network Structure

### Number of Types and Nodes

``` r

# 3 types with 4 nodes each
net <- simulate_htna(
  n_nodes = 4,
  n_types = 3,
  seed = 42
)

dim(net$matrix)
#> [1] 12 12
names(net$node_types)
#> [1] "Metacognitive" "Cognitive"     "Behavioral"
```

### Custom Type Names

``` r

# Educational levels
net <- simulate_htna(
  n_nodes = 5,
  type_names = c("Macro", "Meso", "Micro"),
  seed = 42
)

names(net$node_types)
#> [1] "Macro" "Meso"  "Micro"
```

### Custom Learning Categories

``` r

# Use specific categories for each type
net <- simulate_htna(
  n_nodes = 4,
  n_types = 3,
  type_names = c("SelfReg", "Cognitive", "Social"),
  categories = c("metacognitive", "cognitive", "social"),
  seed = 42
)

# Check nodes in each type
net$node_types
#> $SelfReg
#> [1] "Diagnose" "Regulate" "Plan"     "Judge"   
#> 
#> $Cognitive
#> [1] "Summarize"  "Understand" "Classify"   "Process"   
#> 
#> $Social
#> [1] "Present"  "Respond"  "Teach"    "Question"
```

## Controlling Edge Probabilities

### Within vs. Between Type Connectivity

``` r

# High within-type, low between-type connectivity
clustered <- simulate_htna(
  n_nodes = 4,
  n_types = 3,
  within_prob = 0.6,    # 60% within-type edges
  between_prob = 0.1,   # 10% between-type edges
  seed = 42
)

# More integrated network
integrated <- simulate_htna(
  n_nodes = 4,
  n_types = 3,
  within_prob = 0.3,    # 30% within-type
  between_prob = 0.3,   # 30% between-type (equal)
  seed = 42
)
```

### Analyzing Connectivity Patterns

``` r

# Function to analyze within/between type edge counts
analyze_connectivity <- function(net) {
  mat <- net$matrix
  n_types <- length(net$node_types)
  n_per_type <- net$n_nodes_per_type[1]

  within_count <- 0
  between_count <- 0

  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      type_i <- ceiling(i / n_per_type)
      type_j <- ceiling(j / n_per_type)

      if (mat[i, j] > 0) {
        if (type_i == type_j) {
          within_count <- within_count + 1
        } else {
          between_count <- between_count + 1
        }
      }
    }
  }

  list(
    within = within_count,
    between = between_count,
    ratio = within_count / (within_count + between_count)
  )
}

# Compare networks
analyze_connectivity(clustered)
#> $within
#> [1] 21
#> 
#> $between
#> [1] 8
#> 
#> $ratio
#> [1] 0.7241379
analyze_connectivity(integrated)
#> $within
#> [1] 9
#> 
#> $between
#> [1] 23
#> 
#> $ratio
#> [1] 0.28125
```

## Using with tna Package

### Visualization with plot_htna

``` r

library(tna)

net <- simulate_htna(seed = 42)

# Polygon layout (types arranged in circle)
plot_htna(
  weights = net$matrix,
  node_groups = net$node_types,
  layout = "polygon"
)
```

### Visualization with plot_mlna

``` r

# Layer layout (types as horizontal layers)
plot_mlna(
  weights = net$matrix,
  layers = net$node_types
)
```

### Fitting Models to Multi-Type Data

``` r

# Generate sequences from HTNA matrix
net <- simulate_htna(seed = 42)

# Create sequences using the transition matrix
sequences <- simulate_sequences(
  trans_matrix = net$matrix,
  init_probs = rep(1/nrow(net$matrix), nrow(net$matrix)),
  n_sequences = 200,
  seq_length = 30
)

# Fit TNA model
model <- fit_network_model(sequences, "tna")

# The model can be visualized with type information
plot_htna(
  weights = extract_transition_matrix(model),
  node_groups = net$node_types
)
```

## Advanced Applications

### Group TNA Networks

For analyzing multiple groups with the same structure:

``` r

# Generate group TNA networks
group_networks <- simulate_group_tna_networks(
  n_groups = 5,
  n_actors = 10,
  n_states = 6,
  use_learning_states = TRUE,
  categories = c("metacognitive", "cognitive"),
  seed = 42
)
```

### Using simulate_tna_matrix()

An alternative function for generating HTNA/MLNA matrices:

``` r

# Generate matrix with node types
mat_result <- simulate_tna_matrix(
  n_states = 15,
  matrix_type = "htna",
  n_types = 3,
  type_names = c("Planning", "Execution", "Evaluation"),
  seed = 42
)
```

## Educational Research Examples

### Example 1: Self-Regulated Learning Study

``` r

# Three-phase SRL model
srl_net <- simulate_htna(
  n_nodes = 4,
  type_names = c("Forethought", "Performance", "Reflection"),
  categories = c("metacognitive", "behavioral", "metacognitive"),
  within_prob = 0.5,
  between_prob = 0.2,
  seed = 42
)

names(srl_net$node_types)
#> [1] "Forethought" "Performance" "Reflection"
srl_net$node_types
#> $Forethought
#> [1] "Diagnose" "Regulate" "Plan"     "Judge"   
#> 
#> $Performance
#> [1] "Review"   "Complete" "Diagram"  "Record"  
#> 
#> $Reflection
#> [1] "Self_assess" "Adapt"       "Reflect"     "Regulate"
```

### Example 2: Collaborative Learning Study

``` r

# Individual vs. group regulation
collab_net <- simulate_htna(
  n_nodes = 5,
  type_names = c("Individual", "Shared", "Social"),
  categories = c("cognitive", "group_regulation", "social"),
  within_prob = 0.4,
  between_prob = 0.25,
  seed = 42
)

collab_net$node_types
#> $Individual
#> [1] "Process"    "Memorize"   "Read"       "Generalize" "Compare"   
#> 
#> $Shared
#> [1] "Coregulate" "Cohesion"   "Plan"       "Adapt"      "Synthesis" 
#> 
#> $Social
#> [1] "Feedback" "Explain"  "Answer"   "Help"     "Critique"
```

### Example 3: Multimodal Data Integration

``` r

# Different data sources
multimodal_net <- simulate_htna(
  n_nodes = 4,
  type_names = c("ClickStream", "Survey", "Observation"),
  categories = c("lms", "affective", "behavioral"),
  seed = 42
)

multimodal_net$node_types
#> $ClickStream
#> [1] "Resource" "Submit"   "View"     "Course"  
#> 
#> $Survey
#> [1] "Manage"     "Interest"   "Embrace"    "Discourage"
#> 
#> $Observation
#> [1] "Record"   "Edit"     "Rehearse" "Write"
```

## Comparing Multi-Type Networks

### Comparing Two HTNA Networks

``` r

# Generate two networks with same structure
net1 <- simulate_htna(n_nodes = 4, n_types = 3, seed = 42)
net2 <- simulate_htna(n_nodes = 4, n_types = 3, seed = 123)

# Compare matrices
cor(as.vector(net1$matrix), as.vector(net2$matrix))
```

### Analyzing Type-Level Patterns

``` r

# Extract submatrices for each type pair
extract_type_block <- function(net, type1, type2) {
  nodes1 <- net$node_types[[type1]]
  nodes2 <- net$node_types[[type2]]
  net$matrix[nodes1, nodes2]
}

net <- simulate_htna(seed = 42)

# Within Metacognitive
meta_block <- extract_type_block(net, "Metacognitive", "Metacognitive")

# Metacognitive to Cognitive
meta_to_cog <- extract_type_block(net, "Metacognitive", "Cognitive")
```

## Summary

### Key Functions

| Function | Purpose |
|----|----|
| [`simulate_htna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md) | Generate multi-type transition matrix |
| [`simulate_mlna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md) | Alias for simulate_htna |
| [`simulate_mtna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md) | Alias for simulate_htna |
| [`simulate_tna_matrix()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_matrix.md) | Alternative matrix generator |

### Key Parameters

| Parameter      | Description                             |
|----------------|-----------------------------------------|
| `n_nodes`      | Nodes per type                          |
| `n_types`      | Number of types                         |
| `type_names`   | Custom type labels                      |
| `within_prob`  | Edge probability within types           |
| `between_prob` | Edge probability between types          |
| `categories`   | Learning state categories for each type |

### Output Components

| Component          | Description                 |
|--------------------|-----------------------------|
| `matrix`           | Transition matrix           |
| `node_types`       | List mapping types to nodes |
| `type_names`       | Vector of type names        |
| `n_nodes_per_type` | Node counts per type        |

## Next Steps

- Explore
  [`vignette("tna-workflow")`](https://pak.dynasite.org/Saqrlab/articles/tna-workflow.md)
  for complete analysis workflows
- See
  [`vignette("bootstrap-power")`](https://pak.dynasite.org/Saqrlab/articles/bootstrap-power.md)
  for power analysis with multi-type networks
- Check `?plot_htna` and `?plot_mlna` in the tna package for
  visualization options
