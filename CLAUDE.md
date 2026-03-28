# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

Saqrlab is an R package (v0.3.0) for simulating Temporal Network Analysis (TNA) data. It generates synthetic datasets with known ground-truth parameters for testing and validating statistical methods. **Network estimation code lives in the sibling package [Nestimate](../Nestimate/)** — Saqrlab is simulation-only.

## Build & Test Commands

```bash
# Run all tests (~275 assertions)
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-simulate_data.R")'

# Rebuild documentation (roxygen2 → man/ + NAMESPACE)
Rscript -e 'devtools::document()'

# Quick package check (skip slow parts)
Rscript -e 'devtools::check(".", args = c("--no-tests", "--no-examples", "--no-vignettes", "--no-manual"))'

# Full package check
Rscript -e 'devtools::check()'

# Install locally
Rscript -e 'devtools::install()'
```

Tests use testthat edition 3. All tests generate mock data internally — no external test data files.

## Architecture

### Simulation Modules (core of the package)

**TNA-specific simulation** — generates fitted `tna` model objects directly:
- `simulate_tna_network.R` — single fitted TNA model with learning states
- `simulate_tna_networks.R` — multiple TNA models, group networks, HTNA/MLNA matrices (`simulate_tna_networks()`, `simulate_group_tna_networks()`, `simulate_htna()`, `simulate_mlna()`, `simulate_mtna()`)
- `simulate_sequences.R` — Markov chain sequences from transition matrices (`simulate_sequences()`)
- `simulate_sequences_advanced.R` — sequences with stability/instability patterns
- `simulate_matrix.R` — raw transition matrices (`simulate_matrix()`, `simulate_tna_matrix()`)

**Generic statistical simulation** — returns `list(data, params)` for ground-truth access:
- `simulate_data.R` — 7 dataset types: ttest, anova, correlation, clusters, factor_analysis, prediction, mlvar. Supports `complexity` parameter for edge-case injection (NA, outliers, ties, etc.) and batch generation.
- `simulate_latent.R` — `simulate_lpa()`, `simulate_lca()`, `simulate_regression()`, `simulate_fa()`, `simulate_seq_clusters()`. All return `list(data, params)`.

**Network/graph simulation**:
- `simulate_edge_list.R`, `simulate_igraph.R`, `simulate_network.R` — social network structures
- `simulate_long_data.R` — hierarchical long-format group data
- `simulate_onehot_data.R` — one-hot encoded sequences

### Comparison & Validation Pipeline

Orchestrates TNA model comparisons using the `tna` package:
- `compare_network_estimation.R` / `compare_estimation.R` — compare TNA model types (tna/ftna/ctna/atna) via sampling stability
- `compare_reliability.R` — reliability analysis across conditions
- `network_comparison.R` — compare two fitted networks (correlation, RMSE, edge recovery)
- `fit_network_model.R` — fit TNA/fTNA/cTNA/aTNA models to sequences
- `run_bootstrap_simulation.R`, `run_network_simulation.R`, `run_grid_simulation.R` — simulation orchestration
- `run_sampling_analysis.R` — sampling distribution and stability analysis

### Data & Reference

- `learning_states.R` — 180+ learning action verbs across 8 categories (metacognitive, cognitive, behavioral, social, motivational, affective, group_regulation, lms). Accessed via `LEARNING_STATES` or `get_learning_states()`.
- `global_names.R` — 300 diverse names for actor simulations. Accessed via `GLOBAL_NAMES` or `get_global_names()`.
- `generate_probabilities.R` — random row-stochastic transition matrices and initial probability vectors.

### S3 Classes

`tna_reliability_comparison`, `network_estimation`, `sim_batch_type`, `sim_batch_full` — each with print/plot methods.

### Data Flow

Simulation (`simulate_*.R`) → TNA model fitting (`fit_network_model()` via `tna` package) → comparison (`compare_*`, `run_*`) → visualization (`plot_tna_comparison()`, `plot_network_estimation()`).

## Key Design Decisions

- **list(data, params) pattern**: All `simulate_latent.R` functions return `list(data, params)` so ground-truth parameters are JSON-serializable and testable.
- **simulate_data() seed-driven structure**: The seed determines both the random data AND structural parameters (n, effect sizes, n_groups), so each seed produces a structurally unique dataset.
- **simulate_data() complexity injection**: Edge cases (NA, outliers, constant columns, etc.) are injected post-generation. In batch mode with `complexity = "auto"`, 0–3 cases are randomly selected per dataset (seed-driven, reproducible).
- **Package split**: Network estimation (build_network, bootstrap_network, permutation_test, GIMME, MCML, HON, temporal_network, etc.) moved to [Nestimate](../Nestimate/). Saqrlab depends on `tna` for model fitting but does no estimation itself.
- **Learning states**: 8 categories with realistic verb names. Used by `simulate_sequences()` and `simulate_tna_network()` to generate human-readable state labels.

## Project Workflow Files

- `LEARNINGS.md` — accumulated implementation insights (read before starting work)
- `HANDOFF.md` — current session state and next steps (read at session start, overwrite at end)
- `CHANGES.md` — chronological changelog, newest first
- `FEATURES.md` — concise feature list (update when adding features)

## Dependencies

Core imports: tna, seqHMM, dplyr, tidyr, igraph, ggplot2, parallel, future, future.apply, progressr, lhs. Suggested: testthat, knitr, rmarkdown, network, sna, ggridges, cograph, glmnet, lavaan.
