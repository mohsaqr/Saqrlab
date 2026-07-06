# Changelog

## Saqrlab 0.4.0

Saqrlab is now **simulation-only**: network *estimation* code
(build_network, bootstrap_network, permutation_test, temporal_network,
mlvar, GIMME, MCML, higher-order networks, etc.) has moved to the
sibling package [Nestimate](https://github.com/mohsaqr/Nestimate). This
release focuses the package on synthetic-data generation with known
ground truth.

### New Features

#### Unified `saqr_sim` interface

- New `saqr_sim` S3 class wraps every numerical simulator’s output with
  a consistent `$data` / `$params` / `$type` / `$seed` interface, plus
  `[`, [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`head()`](https://rdrr.io/r/utils/head.html),
  [`tail()`](https://rdrr.io/r/utils/head.html),
  [`dim()`](https://rdrr.io/r/base/dim.html),
  [`names()`](https://rdrr.io/r/base/names.html),
  [`str()`](https://rdrr.io/r/utils/str.html),
  [`print()`](https://rdrr.io/r/base/print.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods. `$data` /
  `$params` access is backward-compatible.
- [`simulate()`](https://pak.dynasite.org/Saqrlab/reference/simulate.md)
  — single dispatcher that routes a `type` string to the matching
  explicit-parameter simulator.
- [`list_simulators()`](https://pak.dynasite.org/Saqrlab/reference/list_simulators.md)
  — catalogs the available simulation types as a data.frame.
- [`validate_recovery()`](https://pak.dynasite.org/Saqrlab/reference/validate_recovery.md)
  — compares a method’s point estimates to the known true `params`,
  returning a tidy per-parameter recovery table (absolute/relative error
  and a within-tolerance flag).

#### Explicit-parameter simulators

- [`simulate_ttest()`](https://pak.dynasite.org/Saqrlab/reference/simulate_ttest.md),
  [`simulate_anova()`](https://pak.dynasite.org/Saqrlab/reference/simulate_anova.md),
  [`simulate_correlation()`](https://pak.dynasite.org/Saqrlab/reference/simulate_correlation.md),
  [`simulate_clusters()`](https://pak.dynasite.org/Saqrlab/reference/simulate_clusters.md),
  [`simulate_prediction()`](https://pak.dynasite.org/Saqrlab/reference/simulate_prediction.md)
  — generate data from user-specified means, SDs, correlation/covariance
  matrices, cluster centers, and regression coefficients; ground-truth
  effect sizes returned in `params`.
- [`simulate_longitudinal()`](https://pak.dynasite.org/Saqrlab/reference/simulate_longitudinal.md)
  — VAR(1) multilevel panel data for mlVAR/ESM with explicit temporal,
  contemporaneous, and between-person matrices.

#### New simulation modules

- [`inject_missingness()`](https://pak.dynasite.org/Saqrlab/reference/inject_missingness.md)
  — adds MCAR, MAR, or MNAR missingness to any simulated data frame at a
  controlled rate.
- [`simulate_mlm()`](https://pak.dynasite.org/Saqrlab/reference/simulate_mlm.md)
  and
  [`simulate_growth()`](https://pak.dynasite.org/Saqrlab/reference/simulate_growth.md)
  — multilevel and latent-growth / longitudinal data from explicit fixed
  effects and variance components.
- [`simulate_irt()`](https://pak.dynasite.org/Saqrlab/reference/simulate_irt.md)
  — item response theory data (1PL/2PL/3PL) from explicit item
  parameters and latent abilities.
- [`simulate_survival()`](https://pak.dynasite.org/Saqrlab/reference/simulate_survival.md)
  — time-to-event data with explicit hazard, covariate effects, and
  censoring.
- [`simulate_hmm()`](https://pak.dynasite.org/Saqrlab/reference/simulate_hmm.md)
  — hidden Markov model sequences from explicit emission and transition
  probabilities.

#### Scenario presets

- [`list_scenarios()`](https://pak.dynasite.org/Saqrlab/reference/list_scenarios.md),
  [`get_scenario()`](https://pak.dynasite.org/Saqrlab/reference/get_scenario.md),
  and
  [`run_scenario()`](https://pak.dynasite.org/Saqrlab/reference/run_scenario.md)
  — named simulation recipes for common designs.
- [`tidy_simulation_results()`](https://pak.dynasite.org/Saqrlab/reference/tidy_simulation_results.md)
  — flattens simulation results to a data.frame.
- [`export_simulation()`](https://pak.dynasite.org/Saqrlab/reference/export_simulation.md)
  — writes results to disk.

### Testing

- Added a cross-package fixture-contract test that guards the simulation
  output contract relied on by downstream packages.

### Bug Fixes

- Fixed
  [`simulate_group_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_group_tna_networks.md)
  (group TNA network generation).
- Fixed
  [`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md)
  / network-comparison summaries.
- Fixed grid-simulation result summaries.

## Saqrlab 0.1.0

Initial release of Saqrlab - Simulation and Analysis Tools for Temporal
Network Analysis.

### New Features

#### Data Simulation

- [`simulate_matrix()`](https://pak.dynasite.org/Saqrlab/reference/simulate_matrix.md) -
  Generate simple transition matrices with learning state names
- [`simulate_htna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md),
  [`simulate_mlna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md),
  [`simulate_mtna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md) -
  Multi-type matrices for hierarchical/multilevel network analysis
- [`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md) -
  Generate Markov chain sequences with optional learning states
- [`simulate_sequences_advanced()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences_advanced.md) -
  Sequences with stability modes for realistic patterns
- [`simulate_long_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_long_data.md) -
  Hierarchical group data with actors, groups, and courses
- [`simulate_onehot_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_onehot_data.md) -
  One-hot encoded format for sequence data
- [`simulate_edge_list()`](https://pak.dynasite.org/Saqrlab/reference/simulate_edge_list.md) -
  Social network edge lists simulation

#### Network Generation

- [`simulate_tna_datasets()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md) -
  Complete TNA datasets with sequences and parameters
- [`simulate_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md) -
  Fitted TNA models with random parameters
- [`simulate_group_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_group_tna_networks.md) -
  Group TNA models for multi-group analysis
- [`simulate_tna_matrix()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_matrix.md) -
  HTNA/MLNA matrices with node types
- [`generate_probabilities()`](https://pak.dynasite.org/Saqrlab/reference/generate_probabilities.md) -
  Random transition matrices and initial probabilities

#### Model Fitting

- [`fit_network_model()`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md) -
  Fit TNA, fTNA, cTNA, and aTNA models

#### Network Comparison

- [`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md) -
  Compare two networks using correlation, RMSE, MAE, and more
- [`compare_centralities()`](https://pak.dynasite.org/Saqrlab/reference/compare_centralities.md) -
  Compare centrality profiles between networks
- [`compare_edge_recovery()`](https://pak.dynasite.org/Saqrlab/reference/compare_edge_recovery.md) -
  Edge recovery metrics (precision, recall, F1)

#### Batch Processing

- [`batch_fit_models()`](https://pak.dynasite.org/Saqrlab/reference/batch_fit_models.md) -
  Fit models to multiple datasets in parallel
- [`batch_apply()`](https://pak.dynasite.org/Saqrlab/reference/batch_apply.md) -
  Apply any function to a list of objects

#### Bootstrap & Simulation Studies

- [`run_bootstrap_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_simulation.md) -
  Bootstrap analysis for stability testing
- [`run_grid_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_grid_simulation.md) -
  Parameter grid search simulations
- [`run_network_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_network_simulation.md) -
  Comprehensive model comparison studies
- [`run_bootstrap_iteration()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_iteration.md) -
  Evaluate single bootstrap runs
- [`summarize_grid_results()`](https://pak.dynasite.org/Saqrlab/reference/summarize_grid_results.md) -
  Analyze grid simulation output

#### Learning States

- [`get_learning_states()`](https://pak.dynasite.org/Saqrlab/reference/get_learning_states.md) -
  Get learning verbs by category
- [`list_learning_categories()`](https://pak.dynasite.org/Saqrlab/reference/list_learning_categories.md) -
  Show available categories with counts
- [`select_states()`](https://pak.dynasite.org/Saqrlab/reference/select_states.md) -
  Intelligent state selection based on network size
- `LEARNING_STATES` - Full dataset of 180+ learning verbs in 8
  categories
- `GLOBAL_NAMES` - 300 diverse names for simulation
- [`get_global_names()`](https://pak.dynasite.org/Saqrlab/reference/get_global_names.md) -
  Retrieve names from the global names dataset

#### Utilities

- [`generate_param_grid()`](https://pak.dynasite.org/Saqrlab/reference/generate_param_grid.md) -
  Create parameter combination grids
- [`validate_sim_params()`](https://pak.dynasite.org/Saqrlab/reference/validate_sim_params.md) -
  Validate and set defaults for simulation parameters
- [`summarize_simulation()`](https://pak.dynasite.org/Saqrlab/reference/summarize_simulation.md) -
  Summary statistics for simulation results
- [`summarize_networks()`](https://pak.dynasite.org/Saqrlab/reference/summarize_networks.md) -
  Network-level summaries for model lists

### Improvements

- Standardized parameter names across all functions
- `use_learning_states = TRUE` as default for realistic educational
  simulations
- Parallel processing support via `future` and `future.apply`
- Progress reporting with `progressr`

### Dependencies

- Core: `tna`, `seqHMM`
- Data manipulation: `dplyr`, `tidyr`
- Parallel processing: `parallel`, `future`, `future.apply`
- Utilities: `progressr`, `lhs`
