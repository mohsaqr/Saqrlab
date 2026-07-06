# Simulate a Single TNA Network

Generate a single fitted TNA network model with learning state names.
This is the simplest way to get a ready-to-use tna model object.

## Usage

``` r
simulate_tna_network(
  n_states = 9,
  n_sequences = NULL,
  seq_length = 25,
  categories = c("metacognitive", "cognitive"),
  seed = NULL
)
```

## Arguments

- n_states:

  Integer. Number of states (8-10 recommended). Default: 9.

- n_sequences:

  Integer or NULL. Number of sequences to simulate. If NULL (default),
  randomly selects between 600-1200.

- seq_length:

  Integer. Length of each sequence. Default: 25.

- categories:

  Character vector. Learning state categories to use. Options:
  "metacognitive", "cognitive", "behavioral", "social", "motivational",
  "affective", "group_regulation", "lms", or "all". Default:
  c("metacognitive", "cognitive").

- seed:

  Integer or NULL. Random seed for reproducibility. Default: NULL.

## Value

A tna model object (class "tna") ready for use with all tna functions:
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`centralities()`](http://sonsoles.me/tna/reference/centralities.md),
[`communities()`](http://sonsoles.me/tna/reference/communities.md), etc.

## Details

This function provides the simplest interface for generating a TNA
network:

1.  Selects learning state names from specified categories

2.  Generates random transition probabilities

3.  Simulates Markov chain sequences

4.  Fits and returns a tna model

The returned model is 100% compatible with all tna package functions.

## See also

[`simulate_tna_networks`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
for generating multiple networks,
[`simulate_sequences`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
for generating sequences without fitting,
[`fit_network_model`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md)
for fitting models to existing data.

## Examples

``` r
library(Saqrlab)

# Generate a network (simplest usage)
model <- simulate_tna_network(seed = 42)

# Use with tna functions
if (FALSE) { # \dontrun{
library(tna)
plot(model)
centralities(model)
communities(model)
} # }

# Custom configuration
model <- simulate_tna_network(
  n_states = 8,
  n_sequences = 300,
  categories = "group_regulation",
  seed = 123
)

# Access model components
model$weights
#>                 Adapt   Cohesion  Consensus     Discuss     Emotion     Monitor
#> Adapt     0.013669821 0.08622503 0.01997897 0.045215563 0.288117771 0.200841220
#> Cohesion  0.247881356 0.41472458 0.03230932 0.026483051 0.054555085 0.019067797
#> Consensus 0.075544174 0.42381562 0.08450704 0.189500640 0.072983355 0.005121639
#> Discuss   0.267469880 0.17831325 0.17590361 0.079518072 0.002409639 0.024096386
#> Emotion   0.001490313 0.21609538 0.06110283 0.096870343 0.067064083 0.122205663
#> Monitor   0.292857143 0.10000000 0.25000000 0.123809524 0.130952381 0.023809524
#> Plan      0.133050847 0.33813559 0.26271186 0.005084746 0.057627119 0.005932203
#> Synthesis 0.053691275 0.01006711 0.09060403 0.016778523 0.091722595 0.110738255
#>                 Plan   Synthesis
#> Adapt     0.30494217 0.041009464
#> Cohesion  0.20338983 0.001588983
#> Consensus 0.11139565 0.037131882
#> Discuss   0.25301205 0.019277108
#> Emotion   0.17883756 0.256333830
#> Monitor   0.03571429 0.042857143
#> Plan      0.02881356 0.168644068
#> Synthesis 0.20022371 0.426174497
model$initial_probs
#> NULL
```
