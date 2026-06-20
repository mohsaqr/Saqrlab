# Saqrlab 0.4.0

Saqrlab is now **simulation-only**: network *estimation* code (build_network,
bootstrap_network, permutation_test, temporal_network, mlvar, GIMME, MCML,
higher-order networks, etc.) has moved to the sibling package
[Nestimate](https://github.com/mohsaqr/Nestimate). This release focuses the
package on synthetic-data generation with known ground truth.

## New Features

### Unified `saqr_sim` interface
- New `saqr_sim` S3 class wraps every numerical simulator's output with a
  consistent `$data` / `$params` / `$type` / `$seed` interface, plus `[`,
  `as.data.frame()`, `head()`, `tail()`, `dim()`, `names()`, `str()`,
  `print()`, and `summary()` methods. `$data` / `$params` access is
  backward-compatible.
- `simulate()` — single dispatcher that routes a `type` string to the matching
  explicit-parameter simulator.
- `list_simulators()` — catalogs the available simulation types as a data.frame.
- `validate_recovery()` — compares a method's point estimates to the known true
  `params`, returning a tidy per-parameter recovery table (absolute/relative
  error and a within-tolerance flag).

### Explicit-parameter simulators
- `simulate_ttest()`, `simulate_anova()`, `simulate_correlation()`,
  `simulate_clusters()`, `simulate_prediction()` — generate data from
  user-specified means, SDs, correlation/covariance matrices, cluster centers,
  and regression coefficients; ground-truth effect sizes returned in `params`.
- `simulate_longitudinal()` — VAR(1) multilevel panel data for mlVAR/ESM with
  explicit temporal, contemporaneous, and between-person matrices.

### New simulation modules
- `inject_missingness()` — adds MCAR, MAR, or MNAR missingness to any simulated
  data frame at a controlled rate.
- `simulate_mlm()` and `simulate_growth()` — multilevel and latent-growth /
  longitudinal data from explicit fixed effects and variance components.
- `simulate_irt()` — item response theory data (1PL/2PL/3PL) from explicit item
  parameters and latent abilities.
- `simulate_survival()` — time-to-event data with explicit hazard, covariate
  effects, and censoring.
- `simulate_hmm()` — hidden Markov model sequences from explicit emission and
  transition probabilities.

### Scenario presets
- `list_scenarios()`, `get_scenario()`, and `run_scenario()` — named simulation
  recipes for common designs.
- `tidy_simulation_results()` — flattens simulation results to a data.frame.
- `export_simulation()` — writes results to disk.

## Testing
- Added a cross-package fixture-contract test that guards the simulation output
  contract relied on by downstream packages.

## Bug Fixes
- Fixed `simulate_group_tna_networks()` (group TNA network generation).
- Fixed `compare_networks()` / network-comparison summaries.
- Fixed grid-simulation result summaries.


# Saqrlab 0.1.0

Initial release of Saqrlab - Simulation and Analysis Tools for Temporal Network Analysis.

## New Features

### Data Simulation
- `simulate_matrix()` - Generate simple transition matrices with learning state names
- `simulate_htna()`, `simulate_mlna()`, `simulate_mtna()` - Multi-type matrices for hierarchical/multilevel network analysis
- `simulate_sequences()` - Generate Markov chain sequences with optional learning states
- `simulate_sequences_advanced()` - Sequences with stability modes for realistic patterns
- `simulate_long_data()` - Hierarchical group data with actors, groups, and courses
- `simulate_onehot_data()` - One-hot encoded format for sequence data
- `simulate_edge_list()` - Social network edge lists simulation

### Network Generation
- `simulate_tna_datasets()` - Complete TNA datasets with sequences and parameters
- `simulate_tna_networks()` - Fitted TNA models with random parameters
- `simulate_group_tna_networks()` - Group TNA models for multi-group analysis
- `simulate_tna_matrix()` - HTNA/MLNA matrices with node types
- `generate_probabilities()` - Random transition matrices and initial probabilities

### Model Fitting
- `fit_network_model()` - Fit TNA, fTNA, cTNA, and aTNA models

### Network Comparison
- `compare_networks()` - Compare two networks using correlation, RMSE, MAE, and more
- `compare_centralities()` - Compare centrality profiles between networks
- `compare_edge_recovery()` - Edge recovery metrics (precision, recall, F1)

### Batch Processing
- `batch_fit_models()` - Fit models to multiple datasets in parallel
- `batch_apply()` - Apply any function to a list of objects

### Bootstrap & Simulation Studies
- `run_bootstrap_simulation()` - Bootstrap analysis for stability testing
- `run_grid_simulation()` - Parameter grid search simulations
- `run_network_simulation()` - Comprehensive model comparison studies
- `run_bootstrap_iteration()` - Evaluate single bootstrap runs
- `summarize_grid_results()` - Analyze grid simulation output

### Learning States
- `get_learning_states()` - Get learning verbs by category
- `list_learning_categories()` - Show available categories with counts
- `select_states()` - Intelligent state selection based on network size
- `LEARNING_STATES` - Full dataset of 180+ learning verbs in 8 categories
- `GLOBAL_NAMES` - 300 diverse names for simulation
- `get_global_names()` - Retrieve names from the global names dataset

### Utilities
- `generate_param_grid()` - Create parameter combination grids
- `validate_sim_params()` - Validate and set defaults for simulation parameters
- `summarize_simulation()` - Summary statistics for simulation results
- `summarize_networks()` - Network-level summaries for model lists

## Improvements

- Standardized parameter names across all functions
- `use_learning_states = TRUE` as default for realistic educational simulations
- Parallel processing support via `future` and `future.apply`
- Progress reporting with `progressr`

## Dependencies

- Core: `tna`, `seqHMM`
- Data manipulation: `dplyr`, `tidyr`
- Parallel processing: `parallel`, `future`, `future.apply`
- Utilities: `progressr`, `lhs`
