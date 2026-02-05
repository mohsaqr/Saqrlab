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
- `generate_tna_datasets()` - Complete TNA datasets with sequences and parameters
- `generate_tna_networks()` - Fitted TNA models with random parameters
- `generate_group_tna_networks()` - Group TNA models for multi-group analysis
- `generate_tna_matrix()` - HTNA/MLNA matrices with node types
- `generate_probabilities()` - Random transition matrices and initial probabilities

### Model Fitting & Extraction
- `fit_network_model()` - Fit TNA, fTNA, cTNA, and aTNA models
- `extract_transition_matrix()` - Extract raw or normalized transition matrices
- `extract_initial_probs()` - Extract initial state probabilities
- `extract_edges()` - Extract edge lists with filtering options

### Network Comparison
- `compare_networks()` - Compare two networks using correlation, RMSE, MAE, and more
- `compare_centralities()` - Compare centrality profiles between networks
- `calculate_edge_recovery()` - Edge recovery metrics (precision, recall, F1)

### Data Conversion
- `wide_to_long()` - Convert wide-format sequences to long format
- `long_to_wide()` - Convert long-format data to wide format
- `prepare_for_tna()` - Prepare any data format for tna package
- `action_to_onehot()` - Convert action sequences to one-hot encoding

### Batch Processing
- `batch_fit_models()` - Fit models to multiple datasets in parallel
- `batch_apply()` - Apply any function to a list of objects

### Bootstrap & Simulation Studies
- `run_bootstrap_simulation()` - Bootstrap analysis for stability testing
- `run_grid_simulation()` - Parameter grid search simulations
- `run_network_simulation()` - Comprehensive model comparison studies
- `evaluate_bootstrap()` - Evaluate single bootstrap runs
- `analyze_grid_results()` - Analyze grid simulation output

### Learning States
- `get_learning_states()` - Get learning verbs by category
- `list_learning_categories()` - Show available categories with counts
- `smart_select_states()` - Intelligent state selection based on network size
- `LEARNING_STATES` - Full dataset of 180+ learning verbs in 8 categories
- `GLOBAL_NAMES` - 300 diverse names for simulation
- `get_global_names()` - Retrieve names from the global names dataset

### Utilities
- `create_param_grid()` - Create parameter combination grids
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
