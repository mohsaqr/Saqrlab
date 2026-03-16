# Saqrlab Changes

### 2026-03-16 — Add simulate_fa(); refactor all latent functions to list(data, params) return
- R/simulate_latent.R: Added simulate_fa(loadings, phi, psi, n, seed). Generates MVN data from Sigma = L*Phi*L' + diag(psi) via Cholesky. phi validated via isSymmetric() + chol(). psi auto-computed (1 - communalities) with error on over-factored models. sigma_implied returned in params (raw eigenvalue clamp, no diagonal normalisation). Refactored simulate_lpa(), simulate_lca(), simulate_regression() to return list(data, params) instead of data.frame + attr() — makes ground truth directly accessible to JSON/CLI consumers without R-specific knowledge.
- tests/testthat/test-simulate_latent.R: Updated all 45 existing tests to use r$data / r$params. Added 30 new tests for simulate_fa (structure, column names, params accuracy, sigma_implied correctness, sample covariance recovery, seed reproducibility, oblique model, user-supplied psi, single-factor, and 5 validation error cases). All 75 tests pass.
- docs/simulate_data-cookbook.md: Updated Parameter-Recovery Functions section — new return structure, added simulate_fa() with full parameter table and CLI JSON export example.

### 2026-03-16 — Add simulate_lpa(), simulate_lca(), simulate_regression() for parameter recovery
- R/simulate_latent.R: New file with three parameter-recovery simulation functions. `simulate_lpa(means, sds, props, n, seed)` generates continuous LPA data with explicit profile means (n_vars × n_profiles matrix), SDs, and mixing proportions. `simulate_lca(item_probs, class_probs, n, seed)` generates binary LCA data with explicit item response probability matrix (n_items × n_classes). `simulate_regression(coefs, predictor_sds, error_sd, n, seed)` generates regression data where `lm()` recovers true coefficients. All three store ground truth in `attr(d, "true_params")` and validate inputs strictly. All accept `seed` for reproducibility.
- tests/testthat/test-simulate_latent.R: 45 TDD tests covering structure, binary/numeric columns, class label validity, true_params attributes, parameter recovery at large n, seed reproducibility, and input validation (dimension mismatches, out-of-range probabilities, missing predictor_sds names, non-positive error_sd). All 45 pass.
- docs/simulate_data-cookbook.md: Added "Parameter-Recovery Functions" section covering all three new functions with usage examples, parameter tables, and recovery assertion patterns.

### 2026-03-16 — Remove 20 computation files; enhance simulate_data with complexity + batch
- R/: Removed 20 files moved to Nestimate (boot_glasso, bootstrap_network, build_network, data_conversion, estimate_network, estimator_registry, estimators, extraction, frequencies, gimme, hon, honem, hypa, mcml, mlvar, mogen, permutation_test, utils, temporal_network, velocity_tna)
- R/Saqrlab-package.R: Fixed .onLoad — removed .register_builtin_estimators() call (function now in Nestimate)
- R/simulate_data.R: Added `complexity` parameter ("clean"/"auto"/case vector) with 11 edge cases for enterprise-grade stress testing: na, outliers, ties, duplicates, constant_col, all_na_col (post-injection); tiny_n, heavy_tailed, heteroscedastic, extreme_imbalance, multicollinear (generation-time). Added `n_batch` parameter (placed after `...` to avoid R partial matching of `n=`). Added `type = "batch"` wildcard to generate 1000 datasets of all 7 types. Auto-complexity uses Gumbel-max trick for weighted edge-case sampling. Added print.sim_batch_type and print.sim_batch_full S3 methods. Added %||% operator (moved here from utils.R). `.nearest_pd()` retained locally.
- tests/testthat/test-simulate_data.R: Added 50 new tests (sections: complexity validation, auto complexity, individual edge cases × 11, n_batch single type, type=batch full batch, seed consistency, print methods). All 169 tests pass.

### 2026-03-14 — Package split: computation code → Nestimate
- Created Nestimate package at /Users/mohammedsaqr/Documents/Github/Nestimate/
- Moved 22 R source files (estimation, bootstrap, HON, GIMME, MCML, mlVAR, temporal) to Nestimate
- Moved 18 test files to Nestimate with standalone test helpers
- Nestimate: 2137/2138 tests pass independently
- Saqrlab retains: simulation, testing, verification, data resources
- NOT YET DONE: remove duplicates from Saqrlab, wire up dependency

### 2026-03-14 — Redesign proximity timeline plot
- R/temporal_network.R: Rewrote `.plot_proximity_timeline()` — smooth variable-width micro-segments, Okabe-Ito palette, direct endpoint labels, highlight parameter, alpha=0.4 default. Removed `render_edges`, `smooth` params and `.build_edge_segments()`.
- R/temporal_network.R: Rewrote `.compute_proximity_mds()` — strength-based interleaved slot positioning (hub at y=0, rank-spread to ±1, spline interpolation).
- tests/testthat/test-temporal_network.R: Rewrote Section 19 with 15 tests. All pass.

### 2026-03-11 — Unify HON/HONEM output format with MOGen/HYPA conventions
- R/hon.R: Changed `.hon_sequence_to_node()` from pipe notation (`"B|A.X"`) to readable arrow notation (`"X -> A -> B"`). Updated `.hon_graph_to_edgelist()` to produce unified columns: `path`, `from` (context), `to` (next state), `count`, `probability` (was `weight`), `from_order`, `to_order`. Modified `.hon_extract_rules()` and `.honp_extract_rules()` to return `list(rules, count)` for raw count recovery. Updated `.hon_assemble_output()` to build matrix from graph directly (not edge columns). Updated `summary.saqr_hon()` to parse arrow notation. Updated `plot.saqr_hon()` to use `igraph::graph_from_adjacency_matrix()`.
- R/honem.R: No changes needed — automatically inherits clean node names from HON matrix.
- tests/testthat/test-hon.R: Added `.pipe_to_arrow()` helper for pyHON equivalence tests. Updated all pyHON/pyHON+ equivalence tests to convert Python pipe notation to arrow notation before comparing. Updated column references from `weight` to `probability`. Updated node order parsing from pipe to arrow format. 213 tests pass.
- tests/testthat/test-honem.R: 28 tests pass (no changes needed).
- tests/testthat/test-hypa.R: 32 tests pass (no changes needed).
- Total: 2351 tests pass.

### 2026-03-11 — Add MOGen, HONEM, HYPA higher-order network algorithms
- R/mogen.R: New file (~400 lines). Implements Multi-Order Generative Model (Scholtes 2017; Gote & Scholtes 2023). Builds higher-order De Bruijn graphs at orders k=0..K from sequential trajectory data. `.mogen_count_kgrams()` extracts k-tuples and transitions, `.mogen_transition_matrix()` builds row-stochastic matrices, `.mogen_log_likelihood()` computes hierarchical likelihood (order j for step j, capping at model order k), `.mogen_layer_dof()` counts free parameters per layer. `build_mogen()` selects optimal Markov order via AIC, BIC, or sequential likelihood ratio tests. S3 class `saqr_mogen` with print/summary/plot (IC comparison + likelihood plots).
- R/honem.R: New file (~200 lines). Implements HONEM embedding (Saebi et al. 2020). Takes HON adjacency matrix, builds row-stochastic transition matrix, computes neighborhood matrix S = (1/Z)*sum(exp(-k)*D^{k+1}), then truncated SVD for embeddings. Parameter-free — no random walks or skip-gram. S3 class `saqr_honem` with print/summary/plot (2D embedding scatter).
- R/hypa.R: New file (~250 lines). Implements HYPA path anomaly detection (LaRock et al. 2020). Builds k-th order De Bruijn graph, fits propensity matrix Xi via iterative proportional fitting (Sinkhorn-Knopp) to match observed in/out-strengths, computes hypergeometric CDF p-values per edge. Classifies paths as over/under-represented at user alpha. S3 class `saqr_hypa` with print/summary/plot (expected vs observed scatter).
- tests/testthat/test-mogen.R: 47 tests across 8 sections. Includes order-1 recovery test (200 paths from known Markov chain → BIC selects order 1) and second-order detection test.
- tests/testthat/test-honem.R: 28 tests across 6 sections. Includes cluster separation test (within-cluster nodes closer than across in embedding space).
- tests/testthat/test-hypa.R: 32 tests across 6 sections. Includes IPF convergence test, score bounds, alpha sensitivity test.
- Total: 2350 tests pass (was 2243).

### 2026-03-11 — Add BuildHON+ as method = "hon+" in build_hon()
- R/hon.R: Added 7 new internal functions for HON+ pipeline (~290 lines): `.honp_max_divergence()` (MaxDivergence upper-bound for KLD pruning), `.honp_build_order1()` (eager order-1 with StartingPoints index), `.honp_extend_observation()` (lazy higher-order count builder), `.honp_extend_source_fast()` (lazy cache lookup), `.honp_add_to_rules()` (HON+ variant that triggers lazy extension for missing prefixes), `.honp_extend_rule()` (recursive extension with MaxDivergence pre-check), `.honp_extract_rules()` (top-level HON+ pipeline). Added `method` parameter to `build_hon()` with `"hon+"` as default. Updated roxygen docs with new parameter and Saebi et al. 2020 reference.
- tests/testthat/test-hon.R: Extended from 163 to 212 tests. Added Section 7 (HON+ internals, 13 tests), Section 8 (method parameter, 3 tests), Section 9 (HON vs HON+ equivalence, 3 tests), Section 10 (pyHON+ equivalence via reticulate, 4 tests with exact edge/weight match). Updated Section 6 pyHON tests to use `method = "hon"` explicitly.
- Performance on group_regulation (2000 trajectories, 9 states): HON 0.64s vs HON+ 2.8s at max_order=3. HON+ advantage is parameter-free operation, not raw speed at this scale. HON+ with max_order=99 auto-discovers max order 6 in 2.2s.

### 2026-03-10 — Add build_hon() for Higher-Order Networks
- R/hon.R: New file (~820 lines). Implements BuildHON algorithm (Xu, Wickramarathne & Chawla, 2016) for constructing variable-order dependency networks from sequential trajectory data. `build_hon(data, max_order, min_freq, collapse_repeats)` accepts wide data.frame or list of character vectors. Two-stage pipeline: (1) rule extraction via KL-divergence testing with adaptive threshold `order/log2(1+count)`, recursive order extension, source-suffix cache; (2) network wiring with incoming edge rewiring for higher-order nodes and tail rewiring to highest-order targets. Pipe-notation node naming (e.g., "B|A" = at B, came from A). Returns `saqr_hon` S3 class with adjacency matrix, edge list, node names, order stats. S3 methods: print, summary, plot (igraph). No new dependencies. Faithful line-by-line translation of pyHON.
- R/hon.R (bugfix): `.hon_build_distributions()` now zeros out `count` in place (matching pyHON's `BuildDistributions`), so `.hon_kld_threshold()` uses filtered totals. Previously used unfiltered counts, producing lower thresholds and spurious higher-order rules. Found via group_regulation dataset testing.
- tests/testthat/test-hon.R: New file. 163 tests across 6 sections: input validation (8), observation counting (4), distributions (3), KL-divergence (4), end-to-end pipeline (10), pyHON equivalence via reticulate (5 datasets — exact edge/weight match on all). All pass.

### 2026-03-06 — Add build_gimme() for GIMME (Group Iterative Multiple Model Estimation)
- R/gimme.R: New file. Implements GIMME algorithm using lavaan as SEM backend. `build_gimme()` with simple API: data, vars, id, time, ar, standardize, groupcutoff, seed. Group-level search (modification indices across subjects, prop_cutoff), pruning (remove non-significant group paths), individual-level search (Bonferroni-corrected per person) with testWeights eigenvalue stability check and search-prune-resume cycle. Returns `saqr_gimme` S3 class with temporal/contemporaneous matrices (count + average), per-person standardized coefficients, psi, fit indices, group/individual paths, syntax. S3 methods: print, summary, plot (types: temporal, contemporaneous, individual, counts, fit). Uses cograph for mixed-network plotting. Achieves **exact equivalence** with gimme R package: path counts identical, standardized coefficients diff=0 across 100 randomly generated datasets (70 × 3-var, 30 × 4-var). Three key fixes: standardized betas/z-values (matching gimme's `lavInspect(fit,"std")$beta` and `standardizedSolution()$z`), and resume-after-nonconv in individual search.
- tests/testthat/test-gimme.R: New file. 78 tests across 10 sections: input validation, basic construction, AR paths, path counts, fit indices, reproducibility, S3 methods, gimme package equivalence (path counts identical, coefficients identical, data prep identical), 4-variable test, config storage. All pass.
- NAMESPACE: Added export(build_gimme), S3 methods for saqr_gimme.
- DESCRIPTION: Added lavaan to Suggests.

### 2026-03-05 — Rewrite build_mcml() with native estimation + add bootstrap_mcml()
- R/mcml.R: Rewritten `build_mcml()` to use Saqrlab's own `build_network()` pipeline for between/within estimation instead of `cograph::cluster_summary()`. Between: recodes states to cluster labels, collapses consecutive same-cluster repetitions (via `.collapse_consecutive()`), then calls `build_network()` — produces zero diagonal by construction (only actual cluster switches are counted). Within: filters non-member states to NA then calls `build_network()` per cluster. Fallback for pre-computed inputs (matrix, tna, netobject, edge list, cograph_network): aggregates off-diagonal matrix blocks only (diagonal always zero). `cluster_summary` still computed for `plot_mcml()` compatibility only. 6 data input formats, 5 cluster formats. Added `bootstrap_mcml()` using the same recode/collapse/filter approach via `bootstrap_network()`. S3 class `mcml_bootstrap` with print, summary(level=), plot(type=). `$between` and `$within` are now plain matrices (not tna objects).
- tests/testthat/test-mcml.R: Extended from 61 to 101 tests. Updated section 2 (matrix types), section 7 (sequence vs matrix paths + zero diagonal + collapse behavior), added sections 10-12 (multi-format data, clusters, edge list).
- tests/testthat/test-bootstrap_mcml.R: New file. 56 tests across 10 sections. All pass.
- docs/build_mcml-for-cograph.md: Short reference for cograph developer — clarifies cograph is only used for plotting, not core estimation. Documents collapse-consecutive approach.
- docs/mcml-technical-reference.md: Full technical reference for both functions.

### 2026-03-05 — Add build_mcml() for multi-cluster multi-layer networks
- tmp/test_mcml.Rmd: Visual verification report.

### 2026-02-28 — Add proximity timeline plot for temporal_network
- R/temporal_network.R: Added `.plot_proximity_timeline()` — per-vertex time series plot with configurable Y-axis metric (eigenvector, degree, closeness, betweenness, page_rank, hub_score, authority_score, or 1D MDS proximity). Supports vertex_group coloring, custom vertex_color, labels_at specific bins, smooth (LOESS), render_edges (vertical segments between connected vertices), and auto-hidden legend for >15 vertices. Added `.resolve_proximity_metric()` dispatcher, `.compute_proximity_mds()` for geodesic-distance MDS, and `.build_edge_segments()` helper. Updated `plot.temporal_network()` dispatch to include `"proximity"` type.
- tests/testthat/test-temporal_network.R: Added 12 tests (Section 19) covering all metrics, MDS mode, vertex groups/colors, labels, smooth, edge rendering, legend behavior, and edge cases.
- tmp/test_proximity_timeline.Rmd: Visual verification report rendered to HTML with 16 plot variants.
- Tests: 1793 pass, 0 fail.

### 2026-02-27 — Add technical reports for all modules (FA-style markdown)
- docs/ising-TECHNICAL-REPORT.md: Ising model estimator — nodewise L1-logistic regression, EBIC selection, AND/OR symmetrization, IsingFit exact equivalence (diff=0 across 50+ configs), 15 sections.
- docs/core-architecture-TECHNICAL-REPORT.md: Core architecture — estimator registry pattern, build_network() pipeline, 7 built-in estimators, scaling (4 methods), multilevel decomposition, netobject S3 class, predictability, 15 sections.
- docs/bootstrap-network-TECHNICAL-REPORT.md: Universal bootstrap — fast transition path (pre-computed), full association path, stability/threshold inference, percentile CIs, model pruning, 15 sections.
- docs/permutation-test-TECHNICAL-REPORT.md: Permutation test — edge-level comparison, paired/unpaired, 8 p-value adjustments, Cohen's d effect size, fast transition path, 13 sections.
- docs/temporal-network-TECHNICAL-REPORT.md: Temporal network — temporal BFS, reachability, 21 snapshot metrics, 12 plot types, tsna equivalence across 20+ configs (1e-10), 22 sections.
- docs/velocity-tna-TECHNICAL-REPORT.md: Velocity TNA — regression (OLS/beta/logistic), GLLA, finite difference, effect-size metrics, 3 method agreement within 1e-6, 15 sections.
- docs/simulate-data-TECHNICAL-REPORT.md: Simulate data — 7 dataset types, seed-based variation, base R only, ground truth attributes, 17 sections.
- All rendered to HTML with flatly theme, floating TOC, numbered sections.

### 2026-02-27 — Add mlvar and boot_glasso technical reports (FA-style markdown)
- docs/mlvar-TECHNICAL-REPORT.md: Comprehensive technical report covering 7-step pipeline, fixed-effects OLS with df correction, within-centering equivalence to lmer, EBIC-GLASSO for contemporaneous/between networks, mlVAR cross-validation results (temporal B < 1e-10, p-values < 0.01, residual r > 0.999), 15 sections.
- docs/boot-glasso-TECHNICAL-REPORT.md: Comprehensive technical report covering 4-phase bootstrap pipeline, glassopath optimization, case-dropping CS-coefficient, edge/centrality difference tests, predictability from precision matrix, bootnet cross-validation results (edge CI r > 0.99, CS within 0.15, 3.2x speedup), 19 sections.

### 2026-02-27 — Add boot_glasso scientific technical report
- docs/boot-glasso-technical-report.Rmd: New scientific report matching mlvar style. cograph::splot network visualizations with custom node colors, pie-chart predictability, circular layout, side-by-side original vs thresholded. Dark theme variant. All bootstrap inference plots (edge CIs, inclusion, stability, diff tests). Bootnet equivalence with scatter plots and correlation table. Speed comparison. Code hidden, narrative-focused.
- docs/boot-glasso-technical-report.html: Rendered output.

### 2026-02-27 — Improve boot_glasso diff plots + AKA bootnet equivalence report
- R/boot_glasso.R: Rewrote `.bg_plot_edge_diff()` and `.bg_plot_centrality_diff()` to match bootnet's visual style. Full matrix (both triangles + diagonal) instead of upper-triangle only. Discrete 3-color fill (black = significant, light gray = non-significant, white = diagonal) instead of continuous red→white p-value gradient. Added `order` parameter: `"sample"` (default, sorted by value ascending) or `"id"` (alphabetical). Clean theme with white grid lines between tiles.
- tests/testthat/test-boot_glasso.R: Added 16 new tests for diff plot order parameter, discrete fill, and full matrix dimensions. Total: 244 tests, all pass.
- tmp/test_boot_glasso_aka.Rmd: Added `## Bootnet Equivalence` section comparing boot_glasso vs bootnet on AKA 8-variable real data: original network match, edge CI correlation, strength CI correlation, CS-coefficient comparison, inclusion probability correlation, difference test agreement, speed comparison.

### 2026-02-27 — Add boot_glasso() with bootnet equivalence
- R/boot_glasso.R: New file (~650 lines). Exported `boot_glasso(x, iter, cs_iter, cs_drop, alpha, gamma, nlambda, centrality, cor_method, ncores, seed)` — single-call bootstrap for EBICglasso partial correlation networks. Combines nonparametric edge/centrality CIs, case-dropping CS-coefficient, edge inclusion probabilities, thresholded network, edge/centrality difference tests, and predictability CIs. Accepts data frame, matrix, or glasso netobject. Ten internal helpers. S3 class `boot_glasso` with `print`, `summary(type=)`, and `plot(type=)` methods. 6 plot types: edges, stability, edge_diff, centrality_diff, inclusion, network. Reuses existing GLASSO pipeline. Supports parallel via `mclapply`. No new dependencies.
- R/boot_glasso.R (bootnet equivalence): Rewrote `.bg_case_drop()` to match bootnet's random-per-iteration approach (each iteration randomly picks ONE drop proportion, not fixed iterations per proportion). Changed `.bg_cs_coefficient()` to use `mean(cors > 0.7) > 0.95` (matching bootnet's proportion-above-threshold, not quantile-based). Changed case-dropping correlation from Spearman to Pearson to match bootnet's `cor0()`.
- tests/testthat/test-boot_glasso.R: 228 tests total. Original 217 tests plus 11 new bootnet equivalence tests: original network matches bootnet across 5 seeds (max diff < 0.01), edge CIs r > 0.99, strength CIs r > 0.99, CS-coefficient matches corStability within 0.15, strength matches qgraph::centrality()$InDegree.
- Tests: 1750 full suite, 0 failures, 0 warnings.

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
