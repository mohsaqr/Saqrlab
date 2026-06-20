# Saqrlab — Feature List

> Saqrlab is **simulation-only**. Network *estimation* (build_network,
> bootstrap_network, permutation_test, temporal_network, mlvar, GIMME, MCML,
> higher-order networks, etc.) lives in the sibling package
> [Nestimate](https://github.com/mohsaqr/Nestimate).

## Unified Simulation Interface
- **saqr_sim S3 class** — all numerical simulation functions return `saqr_sim` with `$data`, `$params`, `$type`, `$seed`. Backward-compatible: `$data`/`$params` access unchanged. Supports `[`, `as.data.frame()`, `head()`, `dim()`, `print()`, `summary()`.
- **simulate()** — single dispatcher that routes a `type` string to the matching explicit-parameter simulator; **list_simulators()** catalogs the available types.
- **validate_recovery()** — compares a method's point estimates to the known true `params`, returning a tidy per-parameter recovery table (abs/rel error, within-tolerance flag).
- **Scenario presets** — **list_scenarios()** / **get_scenario()** / **run_scenario()** named simulation recipes; **tidy_simulation_results()** flattens results to a data.frame; **export_simulation()** writes them to disk.

## Data Simulation — Explicit Parameters (saqr_sim return)
- **simulate_ttest()** — two-group data with specified means/SDs; params include Cohen's d
- **simulate_anova()** — multi-group data with specified means/SDs; params include population eta-squared; supports unequal group sizes
- **simulate_correlation()** — multivariate normal from explicit correlation/covariance matrix; custom means and variable names
- **simulate_clusters()** — mixture of Gaussians with specified centers/SDs; per-cluster sizes or proportions; params include true cluster labels
- **simulate_prediction()** — regression with continuous + categorical predictors; explicit coefficients and category effects; params include population R²
- **simulate_longitudinal()** — VAR(1) multilevel panel data for mlVAR/ESM; explicit temporal (B), contemporaneous, between-person matrices; ESM day/beep structure; complexity injection
- **simulate_lpa()** — LPA data with explicit profile means/SDs
- **simulate_lca()** — LCA data with explicit item response probabilities
- **simulate_regression()** — regression data with known named coefficients
- **simulate_fa()** — factor analysis data from explicit loadings/phi/psi; sigma_implied in params
- **simulate_seq_clusters()** — Markov chain sequences with K known-cluster transition matrices

## Data Simulation — Additional Models
- **simulate_mlm()** — multilevel (random-intercept/slope) data with explicit fixed effects and variance components
- **simulate_growth()** — latent growth / longitudinal trajectories with explicit intercept/slope means and covariances
- **simulate_irt()** — item response theory data (1PL/2PL/3PL/GRM) from explicit item parameters and latent abilities
- **simulate_survival()** — survival/time-to-event data with explicit hazard, covariate effects, and censoring
- **simulate_hmm()** — hidden Markov model sequences from explicit emission and transition probabilities

## Missing-Data Injection
- **inject_missingness()** — adds MCAR, MAR, or MNAR missingness to any simulated data frame at a controlled rate

## Data Simulation — Random Parameters
- **simulate_data()** — 15 dataset types: ttest, anova, correlation, clusters, factor_analysis, prediction, mlvar, probit, meta, count, reliability, gam, nma, power, cthmm; 11 complexity cases; batch generation. Seed-driven random structure — returns a bare data.frame

## Visualization
- All S3 classes have print/summary/plot methods
- TNA model-comparison plots (`plot_tna_comparison()`) and sampling-distribution plots (`plot_sampling_distribution()`)
- Network-estimation comparison plots (`plot_network_estimation()`)
