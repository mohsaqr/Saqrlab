# Saqrlab Changes

### 2026-02-27 — Add mlvar() multilevel vector autoregression
- R/mlvar.R: New file (~400 lines). Exported `mlvar(data, vars, id, day, beep, lag, standardize, gamma, nlambda)` estimates 3 networks from ESM/EMA panel data: temporal (directed, within-person OLS with df correction), contemporaneous (undirected, EBIC-GLASSO on residuals), between-subjects (undirected, EBIC-GLASSO on person means). Six internal helpers: `.mlvar_prepare_data()`, `.mlvar_build_lag_pairs()`, `.mlvar_within_center()`, `.mlvar_temporal_ols()`, `.mlvar_contemporaneous()`, `.mlvar_between()`. S3 class `mlvar_result` with `print` and `summary` methods. Reuses existing `.compute_lambda_path()`, `.select_ebic()`, `.precision_to_pcor()` from estimators.R. No new dependencies.
- R/simulate_data.R: Added `"mlvar"` type with `.simulate_mlvar()` generator. Produces panel EMA data with known temporal B matrix, contemporaneous noise structure, and person-specific means. Attributes: `true_temporal`, `true_contemporaneous`, `vars`.
- tests/testthat/test-mlvar.R: New file. 121 tests covering input validation (8), lag pair construction (8), centering (5), temporal OLS (10), contemporaneous network (6), between-subjects network (6), S3 methods (4), simulate_data("mlvar") (6), integration recovery (3), mlVAR package equivalence (6: temporal coefficients exact match across 5+20 seeds, p-values <0.01 diff, residual correlations r>0.999, person mean correlations r>0.99, significance agreement >95%). All pass.
- Tests: 1522 full suite, 0 failures, 0 warnings.

### 2026-02-26 — Add simulate_data() for ready-to-use statistical datasets
- R/simulate_data.R: New file. Exported `simulate_data(type, seed)` with 6 types: ttest, anova, correlation, clusters, factor_analysis, prediction. Six internal generators (`.simulate_ttest()`, `.simulate_anova()`, `.simulate_correlation()`, `.simulate_clusters()`, `.simulate_factor_analysis()`, `.simulate_prediction()`) plus `.nearest_pd()` helper for positive-definite matrix projection. Seed varies n (10–500), effect sizes, groups/variables. Base R + stats only, no new dependencies.
- tests/testthat/test-simulate_data.R: New file. 103 tests covering input validation (3), per-type structure (30), analysis correctness (6), and attributes (3). All pass.
- man/simulate_data.Rd: Auto-generated roxygen docs.
- NAMESPACE: Added `export(simulate_data)`.
- Tests: 1385 full suite, 0 failures, 0 warnings.

### 2026-02-24 — Plot temporal paths tree visualization
- R/temporal_network.R: `temporal_paths()` now returns S3 class `"temporal_paths"` (inherits `data.frame`) with `source` and `start_time` attributes. Added `plot.temporal_paths()` S3 method — igraph tree layout matching tsna's `plot.tPath` style: thick red directed edges with arrival-time labels, source vertex highlighted, hierarchical tree from root. Uses `igraph::layout_as_tree()` + `igraph::plot.igraph()`. All defaults overridable via `...`.
- tests/testthat/test-temporal_network.R: Added section 18 with 10 tests — class/attributes (3), plot for chain/random/star/disconnected/undirected/single-source/edge cases (7). All pass.
- Total: 452 temporal_network tests, 0 failures, 0 warnings.

### 2026-02-23 — Clean and anonymize VerVerVerVer.csv for TNA
- VerVerVerVer.csv: Renamed columns to TNA conventions (name→Code, value→Order, GroupID→Group, author_id_student→Actor, session_nr→Session, sequence→Sequence, sequence_length→Sequence_Length, Ngroup→Group_Size, Time_gap→Time_Gap, new_session→New_Session). Dropped 4 identifying columns (message_id, messageID, session_id, user). Added User column (Course_Actor composite).
- VerVerVerVer.csv: Anonymized all identifiers — Course (3 random plausible names), Group (24 random hex IDs), Actor (76 random hex IDs), Time (shifted by random offset preserving all relative gaps). Code values kept as-is (generic regulatory labels).
- Verification: All 13 checks pass — row count (10543), distributions, group sizes, session structure, time gaps, and no original identifiers in output.

### 2026-02-23 — Snapshot-based temporal network metrics (with tsna parity)
- R/temporal_network.R: Added `.compute_snapshot_metrics()` — computes 21 igraph metrics (13 graph-level, 8 node-level) across all time bin snapshots in a single pass. Graph-level: density_bins, reciprocity, mutuality, dyad_census, transitivity, centralization (degree/betweenness/closeness), n_components, triad_census, mean_distance, diameter, assortativity. Node-level: closeness_snapshot, betweenness_snapshot, eigenvector, page_rank, hub_score, authority_score, constraint, coreness.
- R/temporal_network.R: Changed `.build_static_snapshots()` to point-in-time edge inclusion (`onset <= bin_start & terminus > bin_start`) + `igraph::simplify()` for multi-edge collapse, matching tsna's `network.collapse()`.
- R/temporal_network.R: Directed transitivity uses matrix formula `sum(A² * A) / sum(A²)` (directed weak transitivity) matching `sna::gtrans(measure="weak")`. Vacuous case (no two-paths) returns 1 matching sna convention. Undirected uses `igraph::transitivity(type="global")` with NaN→1.
- R/temporal_network.R: Updated print (mean transitivity + reciprocity), summary (graph-level time series section), plot (4 new types: reciprocity, centralization, eigenvector, dyad_census).
- tests/testthat/test-temporal_network.R: Added section 17 with 32 tests — dimensions (5), known values (10), range/validity (5), tsna equivalence via `.compare_snapshot_metrics()` helper across 6 datasets: chain, star, dense directed, large 30v/100e directed, large 20v/60e undirected, 50v/200e directed (6), plot types (5), additional (1). Exact numerical match on density, reciprocity, mutuality, triad census, and transitivity.
- Total: 432 temporal_network tests, 1262 full suite, all pass, 0 warnings.

### 2026-02-23 — Large-network tsna equivalence tests + undirected BFS fix
- R/temporal_network.R: Fixed undirected temporal BFS — now adds reverse edges before BFS so traversal works both directions. Affects `.compute_all_temporal_paths()` and `temporal_paths()`.
- tests/testthat/test-temporal_network.R: Added section 16 with 11 large-network tsna equivalence tests: 3x directed (50v/200e, 3 seeds), 2x undirected (40v/150e, 50v/200e), dense (30v/400e), sparse (15v/20e), time_interval=3 (30v/100e), very large (100v/500e), long time range (40v/200e, t_max=200), dense undirected near-complete (10v/40e). All compare reachability, degree, formation/dissolution, edge durations, and temporal paths against tsna.
- Total: 360 temporal_network tests, 1190 full suite, all pass, 0 warnings.

### 2026-02-23 — Achieve tsna numerical equivalence for temporal_network
- R/temporal_network.R: 6 fixes for tsna parity:
  1. Reachability now includes self (+1) to match `tsna::tReach()`.
  2. Time bins extended one past t_max so dissolution at terminus is captured.
  3. Dissolution interval changed from (bin_start, bin_end] to [bin_start, bin_end) to match `tsna::tEdgeDissolution()`.
  4. Degree computation changed from interval-overlap to point-in-time at bin_start to match `tsna::tDegree()`.
  5. Vertex name sorting is now numeric-aware (sorts "1","2",..."10" numerically, not alphabetically).
  6. Undirected BFS adds reverse edges so traversal works in both directions.
- tests/testthat/test-temporal_network.R: Added section 15 with 20 tsna equivalence tests across 3 datasets (chain, star, random). Updated 5 existing reachability tests for +1 self inclusion.
- Full suite: 1190 tests, 0 failures, 0 warnings.

### 2026-02-23 — Add temporal network analysis
- R/temporal_network.R: New file (~600 lines). Implements `temporal_network()` constructor with: input validation, vertex canonicalization, time-varying degree (vectorized per-bin), temporal BFS (earliest-arrival paths), forward/backward reachability, temporal closeness/betweenness centrality, edge formation/dissolution, edge durations, inter-event times + burstiness coefficient, temporal density, and static igraph snapshots. Exported helpers: `temporal_paths()`, `extract_snapshot()`. S3 methods: `print.temporal_network()`, `summary.temporal_network()`, `plot.temporal_network()` (8 plot types: degree, formation, reachability, centrality, burstiness, duration, iet, snapshot).
- tests/testthat/test-temporal_network.R: New file. 120 tests across 14 categories: input validation (8), construction (7), time-varying degree (7), temporal BFS/paths (10), reachability (7), temporal closeness (6), temporal betweenness (7), formation/dissolution (6), duration/IET/burstiness (7), temporal density (5), S3 methods (7), snapshot extraction (5), edge cases (5), snapshot list (3).
- tmp/test_temporal_network.Rmd: Visual verification report with chain, star, and random network examples.
- No new dependencies. Uses only igraph (already in Imports) and ggplot2 (already in Imports).
- Tests: 950/950 pass (0 failures, 0 warnings).

### 2026-02-22 — Add Ising model network estimator
- R/estimators.R: Appended 5 functions: `.prepare_ising_input()` (binary data validation/cleaning), `.log1pexp()` (numerically stable softplus), `.ising_nodewise_ebic()` (nodewise L1-logistic regression + EBIC selection), `.symmetrize_ising()` (AND/OR symmetrization), `.estimator_ising()` (main entry point).
- R/estimator_registry.R: Registered `"ising"` in `.register_builtin_estimators()`.
- R/estimate_network.R: Added `isingfit = "ising"` alias in `.resolve_method_alias()`.
- R/build_network.R: Added `ising = "Ising Model Network"` to `method_labels` in `print.netobject()`. Added Ising-specific print block showing gamma, rule, and threshold range.
- DESCRIPTION: Added `glmnet` to Suggests.
- tests/testthat/test-estimator_ising.R: New file. 69 tests covering input validation (9), log1pexp (1), structure (7), algorithm correctness (3), parameters (7), symmetrization unit tests (2), integration with build_network (7), bootstrap (1), permutation_test (1).
- tmp/test_ising.Rmd: Visual verification report.
- R/estimators.R: Fixed EBIC formula from `4*gamma*k*log(p)` to `2*gamma*k*log(p-1)` (nodewise EBIC matching IsingFit). Fixed OR symmetrization from "avg nonzero" to `(W+t(W))/2` (matching IsingFit).
- tests/testthat/test-estimator_ising.R: Added 5 IsingFit equivalence tests (AND, OR, 20 random configs). Total: 74 tests, all pass.
- Validated: exact numerical equivalence with IsingFit::IsingFit across 50+ random configurations (adjacency diff = 0, threshold diff = 0).

### 2026-02-22 — Add velocity_tna() with regression, GLLA, and finite difference methods
- R/velocity_tna.R: New file. Main function `velocity_tna()` with three methods: `"regression"` (default, OLS + beta regression), `"glla"`, `"finite_difference"`. Regression provides per-edge slopes, SEs, t-values, p-values, R². Internal helpers: `.velocity_regression()`, `.fit_ols_edge()`, `.fit_beta_edge()`, `.glla_weight_matrix()`, `.velocity_glla()`, `.velocity_finite_diff()`, `.detect_and_extract_matrices()`, `.validate_matrix_list()`, `.matrices_from_raw_data()`. S3 methods: `print.tna_velocity`, `summary.tna_velocity`, `plot.tna_velocity` (network/series/heatmap). Fixed column-major indexing bug in `.velocity_regression()`.
- tests/testthat/test-velocity_tna.R: New file. 101 tests covering input validation, regression (slope recovery, edge_stats, p-values, R², acceleration, constant series), GLLA weight matrix correctness, linear/quadratic trend recovery, finite differences (forward/central/backward), input format dispatch (list/tna_windows/data.frame), S3 methods, method agreement, edge cases.
- Tests: all 101 pass. devtools::check: 0 errors, 0 warnings, 3 notes (all pre-existing).
