# 00 — Architecture Overview

This document describes the high-level architecture of the Saqrlab network analysis system. A replicating LLM should read this first to understand the design philosophy, file layout, data flow, and core abstractions before diving into implementation details.

---

## 1. Design Philosophy

The system follows three core principles:

1. **Single entry point**: `build_network()` is the only function users call to estimate a network. There is no separate function per method. The `method` parameter selects the estimator.

2. **Composable params**: All method-specific arguments travel in a `params = list(...)` argument. This list is stored in the returned object and can be replayed verbatim by `bootstrap_network()`, `permutation_test()`, or any downstream function — without knowing the method's internals.

3. **Unified S3 class**: Every network is a `netobject` (or `netobject_ml` for multilevel). All downstream code operates on this single class, not on method-specific return types.

---

## 2. File Map

All core R files live in `R/`. Line counts are as of the current version.

| File | Lines | Purpose |
|---|---|---|
| `R/estimators.R` | 787 | All 6 built-in estimator implementations + shared helpers |
| `R/build_network.R` | 532 | `build_network()`, S3 methods (print, plot, predictability) |
| `R/estimate_network.R` | 211 | Deprecated `estimate_network()` wrapper + shared internal helpers (alias resolution, scaling, edge extraction, multilevel decomposition) |
| `R/estimator_registry.R` | 157 | Global estimator registry (register, get, list, remove) |
| `R/bootstrap_network.R` | 599 | `bootstrap_network()` + fast precompute + statistics + S3 methods |
| `R/permutation_test.R` | 580 | `permutation_test()` + fast paths + glassopath optimization + S3 methods |
| `R/frequencies.R` | 363 | `frequencies()` and `convert_sequence_format()` |
| `R/Saqrlab-package.R` | 18 | `.onLoad()` hook, `globalVariables()` declarations |

Supporting files (not covered in this guide): `R/simulate_*.R`, `R/run_*.R`, `R/compare_*.R`, `R/summary_functions.R`, `R/utils.R`, etc.

---

## 3. Data Flow

```
Raw data (data frame or matrix)
    │
    ▼
build_network(data, method, params, scaling, threshold, level, id_col)
    │
    ├─ .resolve_method_alias(method)        → canonical method name
    ├─ get_estimator(method)                → {fn, description, directed}
    ├─ [if level != NULL] .decompose_multilevel()
    ├─ do.call(estimator$fn, c(list(data=data), params))
    │       → returns list(matrix, nodes, directed, cleaned_data, ...extras)
    ├─ .apply_scaling(matrix, scaling)
    ├─ threshold: abs(x) < threshold → 0
    ├─ .extract_edges_from_matrix()
    │
    ▼
netobject (S3 class)
    │
    ├──► bootstrap_network()   → saqr_bootstrap (S3 class)
    │       ├─ transition: precompute per-sequence → colSums resampling
    │       └─ association: full estimator per iteration
    │
    └──► permutation_test()    → saqr_permutation (S3 class)
            ├─ transition: precompute per-sequence → group label shuffle
            └─ association: optimized closures (cor/solve/glassopath)
```

---

## 4. The 6 Estimation Methods

Three method families, six methods total:

### Transition Methods (directed networks)

| Canonical name | Aliases | Output | Estimator function |
|---|---|---|---|
| `"relative"` | `"tna"`, `"transition"` | Row-normalized transition probabilities | `.estimator_relative()` |
| `"frequency"` | `"ftna"`, `"counts"` | Raw integer transition counts | `.estimator_frequency()` |
| `"co_occurrence"` | `"cna"` | Positional co-occurrence counts (undirected) | `.estimator_co_occurrence()` |

All three call `.count_transitions()` or `.count_cooccurrence_*()` internally. Input is sequence data (wide format: rows = sequences, columns = time points).

### Association Methods (undirected networks)

| Canonical name | Aliases | Output | Estimator function |
|---|---|---|---|
| `"cor"` | `"corr"`, `"correlation"` | Thresholded correlation matrix | `.estimator_cor()` |
| `"pcor"` | `"partial"` | Unregularized partial correlations via `solve()` | `.estimator_pcor()` |
| `"glasso"` | `"ebicglasso"`, `"regularized"` | EBIC-selected graphical lasso partial correlations | `.estimator_glasso()` |

All three call `.prepare_association_input()` to clean/validate input. Input is per-observation frequency data (rows = observations, columns = variables) or a pre-computed correlation/covariance matrix.

---

## 5. Method Alias Resolution

Defined in `.resolve_method_alias()` (file: `R/estimate_network.R`, lines 44-62):

```
ebicglasso  → glasso
regularized → glasso
partial     → pcor
correlation → cor
corr        → cor
transition  → relative
tna         → relative
counts      → frequency
ftna        → frequency
cna         → co_occurrence
```

Unknown names pass through unchanged (supporting custom estimators).

---

## 6. The `netobject` S3 Class

Every `build_network()` call returns a `netobject`. Fields:

### Core fields (always present)

| Field | Type | Description |
|---|---|---|
| `$data` | matrix, data.frame, or NULL | **Cleaned** data from estimator. Association: numeric matrix (NAs, non-numeric, zero-variance removed). Transition: original data.frame. Matrix input: NULL. |
| `$matrix` | numeric matrix | The estimated network weight matrix (n_nodes x n_nodes) |
| `$nodes` | character vector | Node names (= row/colnames of matrix) |
| `$directed` | logical | TRUE for transition methods, FALSE for association |
| `$method` | character | Canonical method name (after alias resolution) |
| `$params` | list | The `params` argument as passed (for replay) |
| `$scaling` | character vector or NULL | Scaling applied |
| `$threshold` | numeric | Threshold applied (0 = none) |
| `$n_nodes` | integer | Number of nodes |
| `$n_edges` | integer | Number of non-zero edges |
| `$edges` | data.frame | Columns: `from`, `to`, `weight`. Directed: all non-diagonal non-zero. Undirected: upper triangle only. |
| `$level` | character or NULL | `"between"`, `"within"`, or NULL |

### Method-specific extras (carried from estimator output)

| Field | Present when | Type |
|---|---|---|
| `$cleaned_data` | always | Same as `$data` (legacy name, carried as extra) |
| `$cor_matrix` | cor, pcor, glasso | Full correlation matrix before thresholding |
| `$precision_matrix` | pcor, glasso | Inverse correlation (precision) matrix |
| `$frequency_matrix` | relative, frequency | Raw integer count matrix |
| `$n` | association methods | Sample size (number of rows) |
| `$p` | association methods | Number of variables |
| `$lambda_selected` | glasso | Selected regularization lambda |
| `$ebic_selected` | glasso | EBIC value at selected lambda |
| `$lambda_path` | glasso | Full lambda path tried |
| `$ebic_path` | glasso | EBIC values for each lambda |
| `$gamma` | glasso | EBIC gamma hyperparameter |

---

## 7. The `netobject_ml` S3 Class

Returned when `level = "both"`. Structure:

| Field | Type | Description |
|---|---|---|
| `$between` | netobject | Between-person network (aggregated to person means) |
| `$within` | netobject | Within-person network (person-mean centered) |
| `$method` | character | The estimation method used |

Created by `build_network()` recursively calling itself with `level = "between"` and `level = "within"`.

---

## 8. Key Design Decisions

### Why `glasso` over `glassoFast`?
`glassoFast` is ~2x faster but is NOT available as a WebR (WASM) binary. `glasso` IS on the WebR repo. The package targets web deployment, so `glasso` was chosen despite the speed trade-off. `glassoFast` also lacks a `penalize.diagonal` parameter.

### Why `data.table` for long-format counting?
Long-format sequence data requires group-by operations (split by person/sequence ID, order by time, extract consecutive pairs). `data.table` handles this with vectorized C-level grouping, avoiding R-level `split()` + `lapply()` overhead.

### Why store `$data` as cleaned, not raw?
Downstream functions (`bootstrap_network()`, `permutation_test()`) need to resample or pool the data. If `$data` were raw, every iteration would re-run cleaning (drop NAs, zero-variance, non-numeric). Storing cleaned data lets downstream code use it directly. For association methods, `$data` is the cleaned numeric matrix; for transition methods, it's the input data frame.

### Why `tabulate()` instead of `table()`?
`tabulate()` operates on integer-encoded vectors and returns a fixed-length vector. It avoids factor creation, drops no levels, and is significantly faster than `table()` for counting known-range integer indices.

### Why `byrow = TRUE` in `matrix()` after `tabulate()`?
The flat index `(from_int - 1) * n_states + to_int` encodes pairs in row-major order (all destinations for source 1, then all for source 2, etc.). R's `matrix()` fills column-major by default. Without `byrow = TRUE`, from/to get silently swapped — a critical bug that produces incorrect networks with no error message.

### Why running accumulators instead of storing all permutation matrices?
Storing all permutation results requires O(iter x n_nodes^2) memory. Running counters (`exceed_counts`, `sum_diffs`, `sum_diffs_sq`) require only O(n_nodes^2). For 1000 iterations with 9 nodes, this is 81K vs 81M numbers.

### Why `glassopath()` in permutation tests?
`glassopath()` solves the entire regularization path in a single Fortran call with internal warm-starting. This is 14.5x faster than calling `glasso()` separately for each lambda value. Lambda must still be re-selected via EBIC each iteration — using a fixed lambda drops agreement from 99% to 93%.

---

## 9. Package Initialization

In `R/Saqrlab-package.R` (lines 14-17):

```r
.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators()
}
```

This calls `.register_builtin_estimators()` from `R/estimator_registry.R` (lines 142-157), which registers all 6 built-in estimators in the global `.estimator_registry` environment. The registry is a `new.env(parent = emptyenv())` — a clean environment with no parent scope pollution.

Global NSE symbols are declared to avoid R CMD check NOTEs:

```r
utils::globalVariables(c(
  "density", "mean_val", ...,
  ".grp_key", ".", ".SD", ":=",
  "cr_lower", "cr_upper", "weight_x", "weight_y", ".seq_grp"
))
```

---

## 10. Dependencies

Core network analysis uses minimal dependencies:

| Package | Used for | Required? |
|---|---|---|
| `glasso` | Graphical lasso estimation | Only for method="glasso" |
| `data.table` | Long-format counting, summary construction | Yes (transition + summary) |
| `stats` | `cor()`, `quantile()`, `sd()`, `p.adjust()`, `var()`, `complete.cases()`, `aggregate()`, `ave()` | Yes (base R) |
| `cograph` | Network plotting via `splot()` | Optional (suggested) |
| `dplyr` + `tidyr` | `convert_sequence_format()` only | Only for format conversion |

Note: `determinant()` is from `base`, NOT `stats`. Do not use `@importFrom stats determinant`.

---

## 11. Related Guides

| Guide | Content |
|---|---|
| [01_estimators.md](01_estimators.md) | Line-by-line algorithm for all 6 estimators |
| [02_build_network.md](02_build_network.md) | The entry point function walkthrough |
| [03_bootstrap.md](03_bootstrap.md) | Bootstrap system with fast precompute |
| [04_permutation_test.md](04_permutation_test.md) | Permutation test with glassopath optimization |
| [05_registry_and_s3.md](05_registry_and_s3.md) | Registry CRUD and all S3 methods |
| [06_testing_and_tips.md](06_testing_and_tips.md) | Testing strategy, pitfalls, dos/don'ts |
