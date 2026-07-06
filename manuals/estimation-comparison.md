# Estimation Comparison Functions

Saqrlab provides two functions for comparing TNA model estimation performance. Each serves a different purpose.

## Overview

| Function | Input | Method | Purpose |
|----------|-------|--------|---------|
| `compare_network_estimation()` | Your real data | Sample vs Remaining | Test model stability |
| `compare_estimation()` | Simulated data | Sample vs Ground Truth | Test estimation accuracy |

---

## `compare_network_estimation()`

Compare how different TNA model types perform on **your actual data** using a statistically correct sampling procedure.

### How It Works

1. Takes your sequence data
2. For each iteration:
   - Splits data into **sample** (30%) and **remaining** (70%)
   - Builds models on both independent subsets
   - Compares them using `tna::compare()`
3. Aggregates metrics across all iterations
4. Ranks models by Pearson correlation

This approach is statistically correct because it compares models built on **independent data subsets**, providing a true measure of estimation stability.

### Usage

```r
library(Saqrlab)
library(tna)
data(group_regulation)

# Compare default models (tna vs ftna)
results <- compare_network_estimation(group_regulation)

# Compare all 4 model types
results <- compare_network_estimation(
  data = group_regulation,
  model_types = c("tna", "atna", "ctna", "ftna"),
  iterations = 100,
  sampling_percent = 0.3,
  seed = 42
)

# View results
print(results)
results$ranking
results$winner

# Plot (default: histogram)
plot(results)
plot(results, metric_name = "Spearman")
plot(results, plot_type = "boxplot")
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `data` | required | Data frame with sequence data |
| `model_types` | `c("tna", "ftna")` | Models to compare: "tna", "atna", "ctna", "ftna" |
| `sampling_percent` | `0.3` | Proportion for sample subset (0-1) |
| `iterations` | `100` | Number of sampling iterations |
| `model_scaling` | `NULL` | Named list of scaling per model |
| `seed` | `NULL` | Random seed for reproducibility |
| `verbose` | `TRUE` | Print progress messages |

### Output

Returns a `network_estimation` object with:

- `$individual` - Raw results from each iteration
- `$aggregated` - Mean, SD, median, min, max by metric and model
- `$ranking` - Models ranked by Pearson correlation (best to worst)
- `$winner` - Best performing model
- `$params` - Parameters used

### Plot Types

```r
plot(results, plot_type = "histogram")   # Default - distribution with mean/median
plot(results, plot_type = "boxplot")     # Box plots by model
plot(results, plot_type = "density")     # Overlapping densities
plot(results, plot_type = "bar")         # Mean with error bars
plot(results, plot_type = "ridgeline")   # Ridge line plots (requires ggridges)
```

---

## `compare_estimation()`

Test how well different TNA models **recover known ground truth** using simulated data.

### How It Works

1. For each simulation:
   - Generates random transition probabilities (ground truth)
   - Simulates sequences from those probabilities
   - Fits all specified model types
   - Compares each fitted model to ground truth using `tna::compare()`
2. Aggregates metrics across all simulations
3. Ranks models by Pearson correlation

This tests **estimation accuracy** - how well each model recovers the true underlying structure.

### Usage

```r
library(Saqrlab)

# Compare tna vs ftna (default)
results <- compare_estimation(n_simulations = 100, seed = 42)

# Compare multiple models
results <- compare_estimation(
  models = c("tna", "ftna", "ctna", "atna"),
  n_simulations = 500,
  n_sequences = 200,
  seq_length = 25,
  n_states = 6,
  na_range = c(0, 5),
  seed = 42
)

# View results
results$comparison   # Side-by-side metrics
results$ranking      # Best to worst
results$winner       # Top performer
results$summary      # Full summary statistics
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `models` | `c("tna", "ftna")` | Model types to compare |
| `n_simulations` | `1000` | Number of simulations to run |
| `n_sequences` | `200` | Sequences per simulation |
| `seq_length` | `25` | Maximum sequence length |
| `n_states` | `6` | Number of states |
| `na_range` | `c(0, 5)` | Range of NAs per sequence (varying lengths) |
| `scaling` | `"minmax"` | Scaling for tna::compare() |
| `seed` | `NULL` | Random seed |
| `verbose` | `TRUE` | Print progress |
| `parallel` | `FALSE` | Use parallel processing |
| `cores` | `detectCores() - 1` | Number of cores if parallel |

### Output

Returns a list with:

- `$comparison` - Side-by-side comparison of key metrics
- `$summary` - Full summary (mean, sd, median) by model and metric
- `$raw_results` - All simulation results
- `$ranking` - Models ranked by Pearson correlation
- `$winner` - Best performing model
- `$params` - Parameters used

---

## When to Use Which?

### Use `compare_network_estimation()` when:

- You have **real data** and want to assess model stability
- You want to know which model gives most **consistent** results on your data
- You're deciding which model type to use for your analysis
- You want **publication-ready plots** of model comparison

### Use `compare_estimation()` when:

- You want to test **theoretical accuracy** of different models
- You need to know which model best **recovers true structure**
- You're doing **methodological research** on TNA models
- You want to test performance under different conditions (sample size, sequence length, etc.)

---

## Example: Complete Workflow

```r
library(Saqrlab)
library(tna)

# 1. First, test theoretical accuracy with simulated data
sim_results <- compare_estimation(
  models = c("tna", "ftna", "ctna", "atna"),
  n_simulations = 200,
  seed = 42
)
cat("Best at recovering ground truth:", sim_results$winner, "\n")

# 2. Then, test stability on your actual data
data(group_regulation)
real_results <- compare_network_estimation(
  group_regulation,
  model_types = c("tna", "ftna", "ctna", "atna"),
  iterations = 100,
  seed = 42
)
cat("Most stable on real data:", real_results$winner, "\n")

# 3. Visualize the real data comparison
plot(real_results)
plot(real_results, plot_type = "boxplot")
```

---

## Available Model Types

| Code | Full Name | Description |
|------|-----------|-------------|
| `tna` | Transition Network Analysis | Relative probabilities (rows sum to 1) |
| `ftna` | Frequency TNA | Raw transition counts |
| `ctna` | Co-occurrence TNA | Co-occurrence based weights |
| `atna` | Attention TNA | Attention-weighted transitions |

---

## Metrics Compared

Both functions use `tna::compare()` which provides:

**Correlations:**
- Pearson, Spearman, Kendall, Distance correlation

**Dissimilarities:**
- Euclidean, Manhattan, Bray-Curtis, Canberra, Frobenius

**Similarities:**
- Cosine, Jaccard

The **Pearson correlation** is used for ranking models (higher = better).
