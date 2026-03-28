# Saqrlab — Feature List

## Network Estimation
- **build_network()** — unified entry point for 7 methods: relative (TNA), frequency, co-occurrence, correlation, partial correlation, EBICglasso, Ising model
- **Estimator registry** — extensible pattern; register custom estimators with `register_estimator()`
- **Multilevel decomposition** — between/within network estimation via `level` parameter
- **Predictability** — node-level R² from precision matrix (Haslbeck & Waldorp 2018)

## Bootstrap & Inference
- **bootstrap_network()** — universal bootstrap for all 7 methods; stability/threshold inference, CIs, p-values, pruned networks
- **permutation_test()** — edge-level network comparison; paired/unpaired, 8 p-value corrections, Cohen's d; exact match with NCT and tna packages
- **boot_glasso()** — single-call bootstrap for EBICglasso networks; edge/centrality CIs, CS-coefficient, difference tests, inclusion probabilities; bootnet-equivalent

## Temporal Network Analysis
- **temporal_network()** — dynamic network analysis from edge lists with onset/terminus; temporal BFS, reachability, closeness/betweenness centrality, formation/dissolution, edge durations, IET, burstiness, 21 snapshot metrics, 12 plot types, proximity timeline; tsna-equivalent
- **velocity_tna()** — edge velocity/acceleration via regression (OLS/beta), GLLA, or finite difference

## Multilevel & Panel Models
- **mlvar()** — multilevel VAR for ESM/EMA data; temporal (fixed-effects OLS), contemporaneous (EBICglasso on residuals), between-subjects networks; mlVAR-equivalent
- **build_gimme()** — GIMME (Group Iterative Multiple Model Estimation); lavaan-backed uSEM with group + individual path search, testWeights stability, standardized output; gimme-equivalent (100/100 datasets, diff=0)

## Multi-Cluster Multi-Layer (MCML)
- **build_mcml()** — multi-cluster multi-layer networks; 6 data formats, 5 cluster formats; native estimation via build_network pipeline
- **bootstrap_mcml()** — bootstrap inference for between/within cluster networks

## Higher-Order Networks
- **build_hon()** — Higher-Order Network construction from sequential trajectory data; two methods: `"hon+"` (default, parameter-free BuildHON+ with lazy observation building and MaxDivergence pruning — Saebi, Xu et al. 2020) and `"hon"` (original BuildHON — Xu et al. 2016); variable-order dependency detection via KL-divergence; pyHON/pyHON+-equivalent
- **build_mogen()** — Multi-Order Generative Model (Scholtes 2017; Gote & Scholtes 2023); higher-order De Bruijn graphs at orders k=0..K, AIC/BIC/LRT-based optimal Markov order selection; hierarchical likelihood decomposition
- **build_honem()** — HONEM Higher-Order Network Embedding (Saebi et al. 2020); parameter-free embeddings via exponentially-decaying matrix powers + truncated SVD; accepts HON output or any square matrix
- **build_hypa()** — HYPA path anomaly detection (LaRock et al. 2020); multi-hypergeometric null model on k-th order De Bruijn graphs; iterative proportional fitting for propensity matrix; per-edge p-values for over/under-represented paths

## Unified Simulation Interface
- **saqr_sim S3 class** — all numerical simulation functions return `saqr_sim` with `$data`, `$params`, `$type`, `$seed`. Backward-compatible: `$data`/`$params` access unchanged. Supports `[`, `as.data.frame()`, `head()`, `dim()`, `print()`, `summary()`.

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

## Data Simulation — Random Parameters
- **simulate_data()** — 7 dataset types: ttest, anova, correlation, clusters, factor_analysis, prediction, mlvar; 11 complexity cases; batch generation. Seed-driven random structure — returns bare data.frame

## Visualization
- All S3 classes have print/summary/plot methods
- cograph integration for network plots (splot, plot_mcml)
- ggplot2-based temporal plots (degree, formation, reachability, centrality, burstiness, proximity timeline)
