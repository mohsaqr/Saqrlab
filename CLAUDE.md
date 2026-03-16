# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

Saqrlab is an R package (v0.3.0) for Temporal Network Analysis (TNA). It provides network estimation (6 methods), bootstrap/permutation testing, advanced algorithms (GIMME, MCML, HON/HON+, HONEM, HYPA, MOGen), temporal network analysis, and simulation tools.

## Build & Test Commands

```r
# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-build_network.R")

# Rebuild documentation (roxygen2 → man/ + NAMESPACE)
devtools::document()

# Full package check
devtools::check()

# Install locally
devtools::install()
```

Tests use testthat edition 3. No external test data files — all tests generate mock data internally.

## Architecture

### Estimation Pipeline

`build_network()` is the **single entry point** for all network estimation. It dispatches to one of 6 estimators (relative, frequency, co_occurrence, cor, pcor, glasso, ising) via an **estimator registry** (`.estimator_registry` environment in `R/estimator_registry.R`). Returns `netobject` S3 class. Supports multilevel decomposition (`level = "between"/"within"/"both"` → `netobject_ml`).

### Bootstrap & Inference

- `bootstrap_network()` — universal bootstrap for all 6 methods. Has a **fast precomputed path** for transitions (pre-computes per-sequence count matrices, single `tabulate()` per iteration).
- `permutation_test()` — edge-level network comparison (paired/unpaired, 8 p-value corrections). Shares precomputed infrastructure with bootstrap.
- `boot_glasso()` — specialized bootstrap for EBICglasso (edge/centrality CIs, CS-coefficient, difference tests).

### Advanced Methods

- **GIMME** (`R/gimme.R`): lavaan-based uSEM with group + individual path search via modification indices. Uses `testWeights` eigenvalue stability check. Exact equivalence with the gimme R package.
- **MCML** (`R/mcml.R`): Multi-cluster multi-layer networks. Uses native `build_network()` pipeline, NOT `cograph::cluster_summary()`. Between matrix includes self-loops (DO NOT collapse consecutive same-cluster states — explicit design decision).
- **HON family**: `R/hon.R` (BuildHON/BuildHON+), `R/mogen.R` (MOGen), `R/honem.R` (HONEM embedding), `R/hypa.R` (HYPA anomaly detection). All use unified edge format: columns `path`, `from`, `to`, `count`, `probability`. Node names use arrow notation (`"A -> B -> C"`).

### Temporal Networks

`R/temporal_network.R` — dynamic network analysis with temporal BFS (multi-pass convergence for non-DAG edges), 21 snapshot metrics, 12 plot types. tsna-equivalent reachability/paths.

### Data Flow

Simulation (`simulate_*.R`) → format conversion (`R/frequencies.R`, `R/data_conversion.R`) → estimation (`build_network()`) → inference (`bootstrap_network()`, `permutation_test()`) → analysis/plotting.

### S3 Classes

`netobject`, `netobject_ml`, `saqr_bootstrap`, `saqr_permutation`, `mcml_network`, `mcml_bootstrap`, `temporal_network`, `boot_glasso`, `saqr_gimme`, `saqr_hon`, `saqr_honem`, `saqr_hypa`, `saqr_mogen` — each with print/summary/plot methods.

## Key Design Decisions

- **MCML between-cluster**: Does NOT collapse consecutive same-cluster states. Between matrix includes self-loops, matching TNA behavior exactly.
- **MCML internals**: `$between` and `$within` are plain matrices, not tna objects. Uses `build_network()` pipeline; `cograph::cluster_summary` only for `plot_mcml()`.
- **HON output**: Unified edge columns (`path`, `from`, `to`) with readable arrow notation for node names.
- **Bootstrap optimization**: Pre-computed per-sequence counts for 2.8x speedup over tna::bootstrap.
- **Estimator extensibility**: New methods added via `register_estimator()` without modifying core code.

## Project Workflow Files

- `LEARNINGS.md` — accumulated implementation insights (read before starting work)
- `HANDOFF.md` — current session state and next steps (read at session start, overwrite at end)
- `CHANGES.md` — chronological changelog, newest first
- `FEATURES.md` — concise feature list (update when adding features)

## Dependencies

Core imports: tna, seqHMM, dplyr, tidyr, igraph, ggplot2, glasso, data.table, parallel, future, future.apply, progressr, lhs. Suggested: testthat, lavaan, cograph, glmnet, network, sna.
