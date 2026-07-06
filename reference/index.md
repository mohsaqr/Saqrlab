# Package index

## Getting started

One front door for every simulator, plus recovery scoring.

- [`simulate()`](https://pak.dynasite.org/Saqrlab/reference/simulate.md)
  : Simulate Data via a Unified Dispatcher
- [`list_simulators()`](https://pak.dynasite.org/Saqrlab/reference/list_simulators.md)
  : List Available Simulators
- [`saqr_sim()`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
  : Create a saqr_sim Object
- [`validate_recovery()`](https://pak.dynasite.org/Saqrlab/reference/validate_recovery.md)
  [`print(`*`<recovery_result>`*`)`](https://pak.dynasite.org/Saqrlab/reference/validate_recovery.md)
  [`summary(`*`<recovery_result>`*`)`](https://pak.dynasite.org/Saqrlab/reference/validate_recovery.md)
  : Score Parameter Recovery Against Simulated Ground Truth
- [`tidy_simulation_results()`](https://pak.dynasite.org/Saqrlab/reference/tidy_simulation_results.md)
  : Flatten Simulation Results into One Tidy Data Frame
- [`export_simulation()`](https://pak.dynasite.org/Saqrlab/reference/export_simulation.md)
  : Export Simulation Results to CSV

## Scenario presets

Named, reproducible simulation recipes.

- [`list_scenarios()`](https://pak.dynasite.org/Saqrlab/reference/list_scenarios.md)
  : List Available Scenario Presets
- [`get_scenario()`](https://pak.dynasite.org/Saqrlab/reference/get_scenario.md)
  : Get a Scenario Recipe
- [`run_scenario()`](https://pak.dynasite.org/Saqrlab/reference/run_scenario.md)
  : Run a Whole Scenario

## Statistical simulators

Classic designs with explicit ground-truth effect sizes.

- [`simulate_ttest()`](https://pak.dynasite.org/Saqrlab/reference/simulate_ttest.md)
  : Simulate Two-Group Comparison Data (t-test)
- [`simulate_anova()`](https://pak.dynasite.org/Saqrlab/reference/simulate_anova.md)
  : Simulate Multi-Group Comparison Data (ANOVA)
- [`simulate_correlation()`](https://pak.dynasite.org/Saqrlab/reference/simulate_correlation.md)
  : Simulate Correlated Multivariate Data
- [`simulate_clusters()`](https://pak.dynasite.org/Saqrlab/reference/simulate_clusters.md)
  : Simulate Cluster Data with Known Centers
- [`simulate_prediction()`](https://pak.dynasite.org/Saqrlab/reference/simulate_prediction.md)
  : Simulate Prediction/Regression Data with Known Coefficients
- [`simulate_regression()`](https://pak.dynasite.org/Saqrlab/reference/simulate_regression.md)
  : Simulate Linear Regression Data with Known Coefficients

## Latent-variable simulators

Mixtures, classes, factors, and sequence clusters.

- [`simulate_lpa()`](https://pak.dynasite.org/Saqrlab/reference/simulate_lpa.md)
  : Simulate Latent Profile Analysis Data
- [`simulate_lca()`](https://pak.dynasite.org/Saqrlab/reference/simulate_lca.md)
  : Simulate Latent Class Analysis Data
- [`simulate_fa()`](https://pak.dynasite.org/Saqrlab/reference/simulate_fa.md)
  : Simulate Factor Analysis Data with Known Parameters
- [`simulate_seq_clusters()`](https://pak.dynasite.org/Saqrlab/reference/simulate_seq_clusters.md)
  : Simulate Sequence Data with Known Cluster Structure

## Longitudinal, multilevel & growth

VAR(1) panels, nested data, and growth curves.

- [`simulate_longitudinal()`](https://pak.dynasite.org/Saqrlab/reference/simulate_longitudinal.md)
  : Simulate Longitudinal Panel Data
- [`simulate_mlm()`](https://pak.dynasite.org/Saqrlab/reference/simulate_mlm.md)
  : Simulate Two-Level Multilevel (Hierarchical) Data
- [`simulate_growth()`](https://pak.dynasite.org/Saqrlab/reference/simulate_growth.md)
  : Simulate Latent Growth-Curve Data

## IRT, survival & hidden Markov

Item response, time-to-event, and latent-state processes.

- [`simulate_irt()`](https://pak.dynasite.org/Saqrlab/reference/simulate_irt.md)
  : Simulate Item Response Theory (IRT) Data
- [`simulate_survival()`](https://pak.dynasite.org/Saqrlab/reference/simulate_survival.md)
  : Simulate Survival Data (Cox Proportional Hazards)
- [`simulate_hmm()`](https://pak.dynasite.org/Saqrlab/reference/simulate_hmm.md)
  : Simulate Hidden Markov Model Sequences

## Missing-data mechanisms

Inject MCAR / MAR / MNAR missingness with a known mechanism.

- [`inject_missingness()`](https://pak.dynasite.org/Saqrlab/reference/inject_missingness.md)
  : Inject Missing Values Under a Known Mechanism (MCAR / MAR / MNAR)

## Random-parameter (stress) simulation

Seed-driven structure for robustness testing.

- [`simulate_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
  : Simulate Ready-to-Use Statistical Datasets

## TNA network simulation

Fitted TNA models, group networks, HTNA/MLNA/MTNA structures.

- [`simulate_tna_network()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_network.md)
  : Simulate a Single TNA Network
- [`simulate_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
  [`generate_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_networks.md)
  : Simulate TNA Network Objects
- [`simulate_group_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_group_tna_networks.md)
  [`generate_group_tna_networks()`](https://pak.dynasite.org/Saqrlab/reference/simulate_group_tna_networks.md)
  : Simulate Group TNA Network Objects
- [`simulate_htna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
  [`simulate_mlna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
  [`simulate_mtna()`](https://pak.dynasite.org/Saqrlab/reference/simulate_htna.md)
  : Simulate HTNA/MLNA/MTNA Matrix with Node Types
- [`simulate_tna_datasets()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
  [`generate_tna_datasets()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
  [`generate_sequence_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
  : Simulate TNA Datasets (Sequences + Models + Probabilities)
- [`simulate_tna_matrix()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_matrix.md)
  [`generate_tna_matrix()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_matrix.md)
  : Simulate TNA Transition Matrix with Node Groupings

## Sequences

Markov-chain sequences and sampling from fitted models.

- [`simulate_sequences()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences.md)
  : Simulate Markov Chain Sequences (Basic)
- [`simulate_sequences_advanced()`](https://pak.dynasite.org/Saqrlab/reference/simulate_sequences_advanced.md)
  : Simulate Markov Chain Sequences (Advanced)
- [`simulate_tna_datasets()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
  [`generate_tna_datasets()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
  [`generate_sequence_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_tna_datasets.md)
  : Simulate TNA Datasets (Sequences + Models + Probabilities)
- [`sample_tna()`](https://pak.dynasite.org/Saqrlab/reference/sample_tna.md)
  : Sample and Re-estimate TNA Model

## Matrices & probabilities

Raw transition matrices and stochastic probability vectors.

- [`simulate_matrix()`](https://pak.dynasite.org/Saqrlab/reference/simulate_matrix.md)
  : Simulate Network Matrix
- [`generate_probabilities()`](https://pak.dynasite.org/Saqrlab/reference/generate_probabilities.md)
  : Generate Transition and Initial Probabilities

## Networks & graphs

Social-network structures and format helpers.

- [`simulate_network()`](https://pak.dynasite.org/Saqrlab/reference/simulate_network.md)
  : Simulate statnet Network Object
- [`simulate_igraph()`](https://pak.dynasite.org/Saqrlab/reference/simulate_igraph.md)
  : Simulate igraph Network Object
- [`simulate_edge_list()`](https://pak.dynasite.org/Saqrlab/reference/simulate_edge_list.md)
  : Simulate Social Network Edge List
- [`simulate_long_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_long_data.md)
  : Simulate Long Format Sequence Data
- [`simulate_onehot_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_onehot_data.md)
  : Simulate One-Hot Encoded Sequence Data

## Model fitting

Fit TNA-family models to sequences.

- [`fit_network_model()`](https://pak.dynasite.org/Saqrlab/reference/fit_network_model.md)
  : Fit a Temporal Network Analysis Model
- [`cross_validate_tna()`](https://pak.dynasite.org/Saqrlab/reference/cross_validate_tna.md)
  : Cross-Validate TNA Model Types

## Comparison & evaluation

Compare fitted networks, centralities, and estimators.

- [`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md)
  : Compare Two TNA Networks
- [`compare_centralities()`](https://pak.dynasite.org/Saqrlab/reference/compare_centralities.md)
  : Compare Centrality Profiles
- [`compare_edge_recovery()`](https://pak.dynasite.org/Saqrlab/reference/compare_edge_recovery.md)
  [`calculate_edge_recovery()`](https://pak.dynasite.org/Saqrlab/reference/compare_edge_recovery.md)
  : Calculate Edge Recovery Metrics
- [`compare_estimation()`](https://pak.dynasite.org/Saqrlab/reference/compare_estimation.md)
  : Compare Model Estimation Across Simulations
- [`compare_network_estimation()`](https://pak.dynasite.org/Saqrlab/reference/compare_network_estimation.md)
  [`compare_tna_models()`](https://pak.dynasite.org/Saqrlab/reference/compare_network_estimation.md)
  : Compare Network Estimation Across TNA Model Types
- [`compare_reliability()`](https://pak.dynasite.org/Saqrlab/reference/compare_reliability.md)
  : Compare Reliability Across Data Conditions
- [`print(`*`<network_estimation>`*`)`](https://pak.dynasite.org/Saqrlab/reference/print.network_estimation.md)
  : Print Network Estimation Results
- [`plot(`*`<network_estimation>`*`)`](https://pak.dynasite.org/Saqrlab/reference/plot.network_estimation.md)
  : Plot Network Estimation Comparison Results
- [`print(`*`<tna_reliability_comparison>`*`)`](https://pak.dynasite.org/Saqrlab/reference/print.tna_reliability_comparison.md)
  : Print Method for TNA Reliability Comparison
- [`plot(`*`<tna_reliability_comparison>`*`)`](https://pak.dynasite.org/Saqrlab/reference/plot.tna_reliability_comparison.md)
  : Plot Method for TNA Reliability Comparison

## Simulation studies & bootstrap

Grid experiments, bootstrap stability, and sampling analysis.

- [`run_bootstrap_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_simulation.md)
  : Run Bootstrap Simulation Across Multiple Runs
- [`run_bootstrap_iteration()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_iteration.md)
  [`evaluate_bootstrap()`](https://pak.dynasite.org/Saqrlab/reference/run_bootstrap_iteration.md)
  : Run Bootstrap Iteration for a Single Simulation Run
- [`run_grid_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_grid_simulation.md)
  : Run Grid Search Over Simulation Parameters
- [`run_network_simulation()`](https://pak.dynasite.org/Saqrlab/reference/run_network_simulation.md)
  : Run Network Analysis Simulations
- [`run_sampling_analysis()`](https://pak.dynasite.org/Saqrlab/reference/run_sampling_analysis.md)
  : Run Sampling Analysis on TNA Models
- [`summarize_grid_results()`](https://pak.dynasite.org/Saqrlab/reference/summarize_grid_results.md)
  [`analyze_grid_results()`](https://pak.dynasite.org/Saqrlab/reference/summarize_grid_results.md)
  : Summarize Grid Simulation Results
- [`generate_param_grid()`](https://pak.dynasite.org/Saqrlab/reference/generate_param_grid.md)
  [`create_param_grid()`](https://pak.dynasite.org/Saqrlab/reference/generate_param_grid.md)
  : Generate Parameter Grid for Simulations
- [`batch_fit_models()`](https://pak.dynasite.org/Saqrlab/reference/batch_fit_models.md)
  : Fit Models to Multiple Datasets
- [`batch_apply()`](https://pak.dynasite.org/Saqrlab/reference/batch_apply.md)
  : Apply Function to Multiple Models or Datasets

## Visualization

Plots for comparisons and sampling distributions.

- [`plot_tna_comparison()`](https://pak.dynasite.org/Saqrlab/reference/plot_tna_comparison.md)
  : Plot TNA Model Comparison Results
- [`plot_network_estimation()`](https://pak.dynasite.org/Saqrlab/reference/plot_network_estimation.md)
  : Plot Network Estimation Results
- [`plot_sampling_distribution()`](https://pak.dynasite.org/Saqrlab/reference/plot_sampling_distribution.md)
  : Plot Sampling Distribution for a Single Metric

## Summaries & validation

Summarize results and validate simulation parameters.

- [`summarize_simulation()`](https://pak.dynasite.org/Saqrlab/reference/summarize_simulation.md)
  : Summarize Simulation Results
- [`summarize_networks()`](https://pak.dynasite.org/Saqrlab/reference/summarize_networks.md)
  : Summarize Multiple Networks
- [`validate_sim_params()`](https://pak.dynasite.org/Saqrlab/reference/validate_sim_params.md)
  : Validate and Standardize Simulation Parameters

## Learning states & names

Reference datasets of learning actions and diverse actor names.

- [`LEARNING_STATES`](https://pak.dynasite.org/Saqrlab/reference/learning_states.md)
  : Learning State Verbs for TNA Simulation
- [`GROUP_REGULATION_ACTIONS`](https://pak.dynasite.org/Saqrlab/reference/GROUP_REGULATION_ACTIONS.md)
  : Group Regulation Actions
- [`get_learning_states()`](https://pak.dynasite.org/Saqrlab/reference/get_learning_states.md)
  : Get Learning State Verbs
- [`list_learning_categories()`](https://pak.dynasite.org/Saqrlab/reference/list_learning_categories.md)
  : Get Learning State Categories Summary
- [`select_states()`](https://pak.dynasite.org/Saqrlab/reference/select_states.md)
  [`smart_select_states()`](https://pak.dynasite.org/Saqrlab/reference/select_states.md)
  : Select States for Networks
- [`GLOBAL_NAMES`](https://pak.dynasite.org/Saqrlab/reference/GLOBAL_NAMES.md)
  : Global Names Dataset
- [`get_global_names()`](https://pak.dynasite.org/Saqrlab/reference/get_global_names.md)
  : Get Global Names
- [`list_name_regions()`](https://pak.dynasite.org/Saqrlab/reference/list_name_regions.md)
  : List Available Name Regions
