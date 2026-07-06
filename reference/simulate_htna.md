# Simulate HTNA/MLNA/MTNA Matrix with Node Types

Generate a transition matrix with multiple node types for hierarchical
(HTNA), multilevel (MLNA), or multi-type (MTNA) network analysis. By
default creates a 25-node matrix (5 nodes x 5 types) using learning
category names.

`simulate_htna()`, `simulate_mlna()`, and `simulate_mtna()` are aliases
that produce identical output - use whichever name fits your analysis
context.

## Usage

``` r
simulate_htna(
  n_nodes = 5,
  n_types = 5,
  type_names = c("Metacognitive", "Cognitive", "Behavioral", "Social", "Motivational"),
  within_prob = 0.4,
  between_prob = 0.15,
  weight_range = c(0, 1),
  allow_self_loops = FALSE,
  categories = c("metacognitive", "cognitive", "behavioral", "social", "motivational"),
  seed = NULL
)

simulate_mlna(
  n_nodes = 5,
  n_types = 5,
  type_names = c("Metacognitive", "Cognitive", "Behavioral", "Social", "Motivational"),
  within_prob = 0.4,
  between_prob = 0.15,
  weight_range = c(0, 1),
  allow_self_loops = FALSE,
  categories = c("metacognitive", "cognitive", "behavioral", "social", "motivational"),
  seed = NULL
)

simulate_mtna(
  n_nodes = 5,
  n_types = 5,
  type_names = c("Metacognitive", "Cognitive", "Behavioral", "Social", "Motivational"),
  within_prob = 0.4,
  between_prob = 0.15,
  weight_range = c(0, 1),
  allow_self_loops = FALSE,
  categories = c("metacognitive", "cognitive", "behavioral", "social", "motivational"),
  seed = NULL
)
```

## Arguments

- n_nodes:

  Integer. Number of nodes per type. Default: 5.

- n_types:

  Integer. Number of node types. Default: 5.

- type_names:

  Character vector. Names for node types. Default uses learning
  categories: "Metacognitive", "Cognitive", "Behavioral", "Social",
  "Motivational".

- within_prob:

  Numeric. Probability of edges within each type. Default: 0.4.

- between_prob:

  Numeric. Probability of edges between types. Default: 0.15.

- weight_range:

  Numeric vector of length 2. Range for edge weights. Default: c(0, 1).

- allow_self_loops:

  Logical. Allow diagonal entries. Default: FALSE.

- categories:

  Character vector. Learning state categories for node names. Can be
  single (same for all types) or one per type. Default: matches
  type_names.

- seed:

  Integer or NULL. Random seed. Default: NULL.

## Value

A list containing:

- matrix:

  Numeric transition matrix (rows sum to 1).

- node_types:

  Named list mapping type names to node names (for plot_htna/plot_mlna).

- type_names:

  Character vector of type names.

- n_nodes_per_type:

  Integer vector of node counts per type.

## See also

[`simulate_matrix`](https://pak.dynasite.org/Saqrlab/reference/simulate_matrix.md),
[`simulate_tna_matrix`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_matrix.md)

## Examples

``` r
# Default: 5 types x 5 nodes = 25 node matrix
net <- simulate_htna(seed = 42)
net$matrix
#>            Diagnose Regulate   Plan  Judge Reflect Understand Classify Process
#> Diagnose     0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Regulate     0.3553   0.0000 0.0000 0.1239  0.0000     0.0000   0.0000  0.0000
#> Plan         0.0000   0.0000 0.0000 0.0000  0.0742     0.2517   0.0000  0.1720
#> Judge        0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   1.0000  0.0000
#> Reflect      0.0000   0.4497 0.0000 0.0472  0.0000     0.0000   0.0000  0.0000
#> Understand   0.0000   0.0000 0.0000 0.0000  0.1433     0.0000   0.0000  0.0000
#> Classify     0.0000   0.0938 0.0386 0.0000  0.5398     0.0000   0.0000  0.0000
#> Process      0.0000   0.0000 0.0000 0.0000  0.0000     0.1804   0.0408  0.0000
#> Encode       0.0000   0.0000 0.1825 0.1233  0.0000     0.0000   0.0000  0.0711
#> Abstract     0.0000   0.0000 0.0000 0.0000  0.1050     0.2681   0.0000  0.1195
#> Write        0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.0000  0.1345
#> Review       0.0000   0.1998 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Outline      0.0000   0.1009 0.0000 0.2225  0.0000     0.0000   0.0174  0.0000
#> Revise       0.0065   0.0000 0.0000 0.0000  0.0000     0.2815   0.0000  0.1108
#> Draft        0.0000   0.0000 0.2634 0.1173  0.0000     0.1367   0.0000  0.0000
#> Help         0.2454   0.0340 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Critique     0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Contribute   0.0000   0.0000 0.0000 0.2836  0.0000     0.4019   0.0000  0.0000
#> Present      0.0000   0.0000 0.0000 0.0000  0.0000     0.0590   0.5055  0.0000
#> Seek_help    0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Overcome     0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Accomplish   0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.4745  0.0000
#> Improve      0.0000   0.3168 0.0000 0.0000  0.0000     0.0000   0.0000  0.0000
#> Create       0.0000   0.0000 0.0613 0.0000  0.0000     0.0000   0.4310  0.0000
#> Strive       0.0000   0.0000 0.0000 0.0000  0.0000     0.0000   0.3066  0.0000
#>            Encode Abstract  Write Review Outline Revise  Draft   Help Critique
#> Diagnose   0.4135   0.1031 0.0000 0.0000  0.0000 0.0000 0.4833 0.0000   0.0000
#> Regulate   0.0000   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0667   0.0000
#> Plan       0.0000   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Judge      0.0000   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Reflect    0.0000   0.0000 0.0000 0.0000  0.4232 0.0000 0.0000 0.0000   0.0000
#> Understand 0.0000   0.0000 0.0000 0.6351  0.0000 0.0000 0.0000 0.0000   0.0000
#> Classify   0.0000   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Process    0.0540   0.0350 0.0000 0.2191  0.0000 0.1912 0.0000 0.0000   0.1234
#> Encode     0.0000   0.0000 0.0000 0.1986  0.0000 0.0000 0.0000 0.1794   0.0000
#> Abstract   0.0278   0.0000 0.0000 0.0000  0.0100 0.0000 0.0000 0.1013   0.0000
#> Write      0.0000   0.0000 0.0000 0.4047  0.0000 0.0000 0.0000 0.0771   0.0000
#> Review     0.8002   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Outline    0.0000   0.0000 0.0877 0.0000  0.0000 0.0000 0.3018 0.0000   0.2697
#> Revise     0.1670   0.0000 0.0000 0.0439  0.2234 0.0000 0.0000 0.0330   0.0000
#> Draft      0.0000   0.0000 0.1319 0.2287  0.0000 0.1220 0.0000 0.0000   0.0000
#> Help       0.2699   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Critique   0.0571   0.0000 0.0000 0.0000  0.3740 0.0000 0.3779 0.0630   0.0000
#> Contribute 0.0000   0.0000 0.0000 0.0000  0.0000 0.0015 0.0000 0.0000   0.3130
#> Present    0.0000   0.0000 0.0000 0.0000  0.0214 0.0000 0.0000 0.0000   0.2854
#> Seek_help  0.0000   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Overcome   0.0000   0.0000 0.0000 0.0000  0.0000 0.1552 0.0000 0.0000   0.0000
#> Accomplish 0.0000   0.0000 0.5255 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Improve    0.0000   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.2644   0.0893
#> Create     0.3477   0.0000 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#> Strive     0.0000   0.1698 0.0000 0.0000  0.0000 0.0000 0.0000 0.0000   0.0000
#>            Contribute Present Seek_help Overcome Accomplish Improve Create
#> Diagnose       0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Regulate       0.0000  0.0000    0.0000   0.1782     0.0000  0.2760 0.0000
#> Plan           0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Judge          0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Reflect        0.0000  0.0000    0.0798   0.0000     0.0000  0.0000 0.0000
#> Understand     0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.2216
#> Classify       0.0000  0.0000    0.0000   0.3278     0.0000  0.0000 0.0000
#> Process        0.0465  0.0000    0.0000   0.0000     0.1098  0.0000 0.0000
#> Encode         0.2022  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Abstract       0.0000  0.0000    0.0000   0.0000     0.0000  0.3028 0.0000
#> Write          0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Review         0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Outline        0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Revise         0.0000  0.0000    0.0000   0.1340     0.0000  0.0000 0.0000
#> Draft          0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Help           0.1709  0.0000    0.2654   0.0000     0.0000  0.0000 0.0144
#> Critique       0.0000  0.1279    0.0000   0.0000     0.0000  0.0000 0.0000
#> Contribute     0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Present        0.0000  0.0000    0.0000   0.0000     0.0000  0.1140 0.0000
#> Seek_help      0.0000  0.0000    0.0000   0.0000     1.0000  0.0000 0.0000
#> Overcome       0.0000  0.0000    0.0000   0.0000     0.5243  0.3205 0.0000
#> Accomplish     0.0000  0.0000    0.0000   0.0000     0.0000  0.0000 0.0000
#> Improve        0.0000  0.0000    0.0000   0.3294     0.0000  0.0000 0.0000
#> Create         0.0000  0.0000    0.0000   0.1599     0.0000  0.0000 0.0000
#> Strive         0.0000  0.0000    0.0000   0.2171     0.3065  0.0000 0.0000
#>            Strive
#> Diagnose   0.0000
#> Regulate   0.0000
#> Plan       0.5021
#> Judge      0.0000
#> Reflect    0.0000
#> Understand 0.0000
#> Classify   0.0000
#> Process    0.0000
#> Encode     0.0428
#> Abstract   0.0655
#> Write      0.3838
#> Review     0.0000
#> Outline    0.0000
#> Revise     0.0000
#> Draft      0.0000
#> Help       0.0000
#> Critique   0.0000
#> Contribute 0.0000
#> Present    0.0147
#> Seek_help  0.0000
#> Overcome   0.0000
#> Accomplish 0.0000
#> Improve    0.0000
#> Create     0.0000
#> Strive     0.0000
net$node_types  # List format for plot_htna/plot_mlna
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
#> 

# Custom number of types
net <- simulate_htna(n_nodes = 4, n_types = 3, seed = 42)

# Custom type names
net <- simulate_mlna(
  n_nodes = 6,
  type_names = c("Macro", "Meso", "Micro"),
  seed = 42
)

# Use with tna package:
# plot_htna(net$matrix, net$node_types, layout = "polygon")
# plot_mlna(net$matrix, layers = net$node_types)
```
