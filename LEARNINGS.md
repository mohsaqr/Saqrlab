# Saqrlab Project Learnings

### 2026-06-20

- \[cross-package fixture contract\]:
  [`simulate_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
  is the SINGLE critical-path function for ~1010 JSON fixtures checked
  into JStats/Carm (`<JStats>/validation/r-reference/*-ref.R` call
  `simulate_data(type, seed=i)` for type in
  {ttest,anova,correlation,clusters,factor_analysis,prediction}, default
  complexity, seeds 1..1000; plus a `seed=SEED*100+s, complexity="auto"`
  batch path). Any drift in its RNG stream or structural-param
  derivation silently breaks 100% of their equivalence tests while
  Saqrlab’s own suite stays green. Guarded now by
  tests/testthat/test-fixture-contract.R. NEVER refactor
  simulate_data()’s internal loops/draw-order.
- \[expect_snapshot_value skips on CRAN\]: testthat’s
  `expect_snapshot_value()` auto-skips under R CMD check / non-NOT_CRAN
  runs (“Reason: On CRAN”) and never writes `_snaps`. Useless as an
  always-run contract guard. Use a committed RDS golden +
  `expect_identical` instead
  (tests/testthat/fixtures/fixture-contract-golden.rds).
- \[8 broken exported functions found by characterization tests\]:
  simulate_group_tna_networks, generate_group_tna_networks,
  compare_networks, compare_centralities, compare_edge_recovery,
  calculate_edge_recovery, summarize_grid_results, analyze_grid_results
  all errored on basic usage because they called helpers that were never
  defined (`long_to_wide`, `extract_transition_matrix`,
  `check_val_in_range`, `safe_bind_rows`). Invisible because they had
  zero tests. Use `codetools::checkUsagePackage(all=TRUE)` to enumerate
  ALL undefined-symbol references at once instead of discovering them
  one cascade at a time.
- \[tna::centralities returns a tibble\]: class
  `tna_centralities,tbl_df,...` with a `state` COLUMN and positional
  rownames – so `cent[, measure]` stays a 1-col tibble (a list) and
  [`sd()`](https://rdrr.io/r/stats/sd.html)/[`cor()`](https://rdrr.io/r/stats/cor.html)
  choke (“‘list’ object cannot be coerced to type ‘double’”). Coerce
  with [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  and key rownames by `state` before extracting measure columns.
- \[head/tail/str need importFrom utils\]: defining S3 methods
  `head.saqr_sim`/`tail.saqr_sim`/`str.saqr_sim` with `@export` adds
  `S3method(...)` but NOT the generic import. Without
  `@importFrom utils head tail str` the namespace fails to load under R
  CMD check (“object ‘head’ not found whilst loading namespace”) – works
  under devtools::load_all() because utils is attached.
  dim/names/print/summary/as.data.frame are base generics and need no
  import.
- \[non-ASCII: strings ok, comments warn\]: under `Encoding: UTF-8`,
  non-ASCII in STRING LITERALS is permitted (e.g. accented names in
  global_names.R) but non-ASCII in COMMENTS/code triggers the R CMD
  check “non-ASCII characters” WARNING. Em-dash/ellipsis/arrows in
  roxygen are the usual culprits.
- \[agents: verify code claims before acting\]: parallel audit agents
  flagged two false-positive “bugs” – GROUP_REGULATION_ACTIONS
  “undefined” (it IS defined at simulate_long_data.R) and
  `paste0("S", 1:n_states)` producing “S1.0” (R’s `1:5` is always
  integer, so it’s “S1”). Always confirm a flagged bug against the
  running package before “fixing” working code.

### 2026-03-28

- \[saqr_sim backward compat\]: Adding an S3 class to
  `list(data, params)` via `class(obj) <- c("saqr_sim", "list")` is
  fully backward compatible — `$data`, `$params`,
  [`names()`](https://rdrr.io/r/base/names.html) all keep working. The
  `[.saqr_sim` method delegates to `$data` so subscripting works too.
- \[saqr_sim vs data.frame returns\]: Cannot wrap
  [`simulate_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md)
  in saqr_sim — 169 existing tests treat the result as a bare data.frame
  (`d$group`, `names(d)`, `is.data.frame(d)`). The `$` on a list goes to
  list elements, not data columns. Would need `$.saqr_sim` fall-through
  to `$data` columns, which creates name collision risk with fields like
  “type”. Two tiers is the right design: saqr_sim for explicit-parameter
  functions, bare data.frame for random-parameter
  [`simulate_data()`](https://pak.dynasite.org/Saqrlab/reference/simulate_data.md).
- \[nrow/ncol not S3 generics\]: In base R,
  [`nrow()`](https://rdrr.io/r/base/nrow.html) and
  [`ncol()`](https://rdrr.io/r/base/nrow.html) call
  [`NROW()`](https://rdrr.io/r/base/nrow.html)/[`NCOL()`](https://rdrr.io/r/base/nrow.html)
  which use [`dim()`](https://rdrr.io/r/base/dim.html). Defining
  `dim.saqr_sim()` is sufficient — explicit `nrow.saqr_sim` and
  `ncol.saqr_sim` get exported as regular functions (not methods) and
  pollute NAMESPACE. Just define `dim`.
- \[population eta-squared\]: True population eta-squared = ss_between /
  (ss_between + ss_within) where ss_between = sum(n_k \* (mu_k -
  grand_mean)^2) and ss_within = sum((n_k - 1) \* sigma_k^2). This is
  the parameter, not the sample estimate.
- \[population R-squared\]: True R² = var(signal) / (var(signal) +
  error_sd^2). Computed from the actual linear combination variance on
  the generated data, not from the population formula. Accounts for both
  continuous and categorical predictor contributions.
- \[pooled SD for Cohen’s d\]: Use
  `sqrt(((n_a-1)*sd_a^2 + (n_b-1)*sd_b^2) / (n_a+n_b-2))` for the pooled
  SD in Cohen’s d. This is the standard pooled estimator, not the simple
  average of variances.
- \[VAR stationarity\]: A VAR(1) process y(t) = mu + B(y(t-1) - mu) + e
  is stationary iff all eigenvalues of B have modulus \< 1. Rescaling
  `B <- B * (0.95 / max(Mod(eigen(B)$values)))` guarantees this while
  preserving relative coefficient structure.
- 
- \[.ensure_pd vs .nearest_pd\]: Two PD-fixing helpers now exist:
  `.ensure_pd()` in simulate_longitudinal.R (eigenvalue clamp only) and
  `.nearest_pd()` in simulate_data.R (eigenvalue clamp + correlation
  normalisation). The latter normalises to a correlation matrix; the
  former preserves scale for covariance matrices.

### 2026-03-17

- \[package split cleanup\]: When splitting packages, check for (1)
  functions called by remaining code that were moved, (2) test files for
  moved functions, (3) unused DESCRIPTION Imports. All three bit us
  after the Saqrlab/Nestimate split.
- 
- \[package structural review\]: Run `grep -rn "function_name" R/`
  across both packages to verify every called function is defined
  somewhere. Don’t rely on git history or memory.

### 2026-03-16

- \[matrix() with list input\]:
  `matrix(list_of_vectors, nrow=..., byrow=TRUE)` creates a list-matrix
  (each cell holds a vector), NOT a numeric matrix. Use
  `matrix(unlist(list_of_vectors), nrow=..., byrow=TRUE)` to get a
  numeric matrix. [`rowSums()`](https://rdrr.io/r/base/colSums.html)
  will error on a list-matrix with “‘x’ must be numeric”.
- \[Reduce accumulate for Markov chains\]:
  `Reduce(f, x=seq_len(T-1), accumulate=TRUE, init=s0)` generates T
  states (init + T-1 more), perfect for Markov chain simulation without
  for loops. Returns a list; use
  [`unlist()`](https://rdrr.io/r/base/unlist.html) if all states are
  scalar strings.
- 

### 2026-02-15

- \[tna input formats\]:
  [`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md)
  accepts both a data frame of sequences (rows=sequences, cols=time
  points) and a square numeric matrix (transition frequencies). Column
  prefix can be “T” (as in `group_regulation`) not just “V”.
- \[group_regulation_long\]: Structure is `Actor` (int), `Achiever`
  (chr), `Group` (num), `Course` (chr), `Time` (POSIXct), `Action`
  (chr). Has 27,533 rows, 9 unique actions (adapt, cohesion, consensus,
  coregulate, discuss, emotion, monitor, plan, synthesis).
- \[group_regulation\]: Wide format with T1-T26 columns (not V1-V26).
  2000 rows. NAs are trailing (variable-length sequences).
- \[tna model structure\]: `tna` object is a list with `weights`
  (transition matrix), `inits` (initial probs), `labels` (state names),
  `data` (encoded sequences). `attr(model, "type")` is “relative” for
  probability weights.
- \[transition counting\]: To build a frequency matrix from sequences,
  count consecutive pairs (action\[t\], action\[t+1\]) within each
  sequence. For long data, split by ID and order by time first. For wide
  data, iterate across columns per row. Use
  [`table()`](https://rdrr.io/r/base/table.html) with factors for
  vectorized counting.
- \[frequencies.R\]: File contains two functions: `frequencies()` builds
  a transition frequency matrix (state x state),
  `convert_sequence_format()` converts to frequency counts, onehot,
  edgelist, or follows format. Both accept wide and long input.
- 
- \[pivot_longer default\]:
  [`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)
  creates a “name” column by default when `names_to` is not specified.
  This column is harmless if not used in downstream operations.
- \[R CMD check\]: Package has pre-existing warning about `\usage`
  entries in some man pages (not from new code). Hidden files (.claude,
  .github) and non-portable file paths in docs/ cache trigger NOTEs.
- 
- \[glasso\]: `glasso::glasso()` returns list with `w` (estimated
  covariance), `wi` (estimated precision/inverse covariance). Warm
  starts use `start = "warm"` with `w.init` and `wi.init`. Cold start is
  default.
- \[EBICglasso validation\]: Our native EBIC + glasso implementation
  matches `qgraph::EBICglasso` exactly on test data (same partial
  correlations, same edges, same sparsity).
- \[compositional frequency data\]: Per-sequence action frequencies are
  compositional (rows sum to ~constant sequence length). This causes all
  partial correlations to be strongly negative. Expected behavior, not a
  bug.
- \[glasso vs glassoFast\]: Both have zero deps and use Fortran.
  `glassoFast` is ~2x faster but NOT available as a WebR (WASM) binary.
  `glasso` IS on the WebR repo. `glassoFast` also lacks
  `penalize.diagonal` param (always penalizes diagonal; workaround: pass
  rho as a matrix with zero diagonal). Chose `glasso` for WebR
  compatibility.
- \[engagement stslist\]:
  [`tna::engagement`](http://sonsoles.me/tna/reference/engagement.md) is
  an stslist object. Converting to frequency via
  `convert_sequence_format(as.data.frame(engagement), format = "frequency")`
  produces a `%` column (void/padding marker). Non-syntactic column
  names like `%` need auto-cleaning.
- \[cograph::splot\]: Preferred over `tplot()` for richer rendering.
  Supports `pie_values`/`pie_colors` for predictability pies,
  `donut_values`/`donut_color` for donuts, `node_fill` for node colors,
  `edge_labels = TRUE` for edge weights, `title` param directly (unlike
  tplot). Use `do.call(cograph::splot, dots)` to pass args
  programmatically. Default layout is `"oval"`.
- \[cograph::tplot\]: Thin qgraph-compatible wrapper around splot. Uses
  `pie`/`pieColor` (qgraph names). Does NOT accept `main` — use `title`
  param or [`title()`](https://rdrr.io/r/graphics/title.html) after
  call.
- \[build_network\]: Renamed from `pcor_network()` to `build_network()`
  with `method` parameter. Supports `"glasso"` (aliases: `"ebicglasso"`,
  `"regularized"`), `"pcor"` (alias: `"partial"`), `"cor"` (alias:
  `"correlation"`). S3 class is `"psych_network"`. Main result field is
  `$network_matrix` (was `$pcor_matrix`). Auto-cleans non-syntactic
  column names, all-NA columns, NA rows, and zero-variance columns. Uses
  `glasso` package (WebR-compatible) for the glasso method.
- \[build_network multilevel\]: `level` parameter (`"between"`,
  `"within"`, `"both"`) decomposes data before estimation. Between =
  aggregate to person means; Within = person-mean center each variable
  (drops persons with \<2 obs). `level = "both"` returns
  `psych_network_ml` S3 class. `id_col` is required when `level` is set.
  Must be data frame input (matrix not supported). The `id_vals`
  alignment with cleaned rows uses a `row_mask` boolean vector.
- 
- 
- \[predictability\]: Node predictability (R²) computed as
  `1 - 1/Omega_jj` for glasso/pcor (from precision matrix). For cor
  method, uses multiple R² from correlation matrix subsets
  (`crossprod(r, solve(R_nn, r))`). Follows Haslbeck & Waldorp (2018).
  Validated against mgm: r=0.999, mean \|diff\|=0.008. Zero dependencies
  added. Shown on plots via `pie_values` in splot, not printed by
  default.
- \[plot styling\]: Use
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  with `pie_values` for predictability, `node_fill` for pastel colors,
  `edge_labels = TRUE`, `theme = "colorblind"`, `node_size = 8`.
  Hand-picked palette of 15 distinct pastels cycles via
  [`rep_len()`](https://rdrr.io/r/base/rep.html). For multilevel, use
  `par(mfrow = c(1, 2))` with
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html) to restore.
- \[group_regulation_long nesting\]: Each Actor appears in exactly 1
  Group (2000 actors, 200 groups, 10 actors/group). For multilevel
  decomposition, use Group as grouping variable
  (`id_col = c("Group", "Actor")` — first element is the grouping
  variable).

### 2026-03-11

- \[HON+ lazy observation\]: BuildHON+ (Saebi et al. 2020) only builds
  order-1 counts eagerly, then extends lazily via StartingPoints index.
  Each source key maps to a list of `c(traj_index, position)` pairs.
  Extension looks one step back in each trajectory to build higher-order
  counts on demand. This avoids building unused high-order observations.
- 
- 
- \[HON+ R performance\]: HON+ is slower than HON in R at moderate scale
  (2.8s vs 0.6s on 2000 trajectories, max_order=3) due to higher
  per-operation overhead of environment-based hash maps and list
  appending for StartingPoints. The advantage is parameter-free
  operation (no max_order tuning needed). At max_order=99
  (parameter-free mode), HON+ auto-discovers the right order (6) in
  2.2s.
- \[reticulate testthat\]: Three reticulate quirks in test_file
  context: (1) `py$x <- val` fails — use
  `py_set_attr(py, "x", val)`, (2) Python list elements don’t
  auto-convert — use `py_to_r()` first, (3) Python tuples become R
  environments — return dicts instead for reliable conversion.
- \[HON+ equivalence threshold\]: HON+ uses strict `>` comparison for
  KLD vs threshold (matching pyHON+). With exactly 3 observations, KLD
  can equal threshold exactly (e.g., both = 1.0), so extension doesn’t
  trigger. Need 4+ observations to push threshold below KLD.
- \[MOGen De Bruijn k-grams\]: At order k, nodes are k-tuples encoded
  with `.HON_SEP`. Edges connect overlapping k-tuples extracted from
  (k+1)-length sub-paths. Use `\x02` as edge pair separator (different
  from `\x01` used for k-gram encoding).
  [`table()`](https://rdrr.io/r/base/table.html) on paste’d keys is the
  most efficient counting approach.
- \[MOGen hierarchical likelihood\]: For order-k model, step j uses
  order min(j-1, k). Step 1 = order 0 (marginal), step 2 = order 1, …,
  step k+1 onward = order k. Cumulative DOF = sum of per-layer DOFs.
  This makes nested model comparison valid for LRT.
- 
- \[HONEM not random walk\]: HONEM (Saebi et al. 2020) does NOT use
  random walks or skip-gram/Word2Vec. It’s direct matrix factorization:
  neighborhood matrix S = (1/Z)*sum(exp(-k)*D^{k+1}) then truncated SVD.
  Parameter-free.
- \[HYPA IPF convergence\]: Iterative proportional fitting
  (Sinkhorn-Knopp) for HYPA’s Xi matrix converges quickly (~10-20
  iterations for tol=1e-6). Must respect edge structure: Xi\[i,j\]=0
  wherever adj\[i,j\]=0. Row/col scaling preserves the sparsity pattern.
- \[HYPA hypergeometric\]: `phyper(k, K, N-K, n)` computes CDF P(X\<=k).
  Must [`round()`](https://rdrr.io/r/base/Round.html) Xi values to
  integers for hypergeometric parameters. Clamping K to \[0, N\]
  prevents invalid parameters.

### 2026-02-16

- \[Bayesian transition networks\]: Not viable without an informative
  prior network. Uniform Dirichlet prior just does Laplace smoothing,
  making all transitions non-zero → bootstrap flags everything as
  significant. Only useful when researcher has a theoretical/prior
  transition matrix to update with data.
  `markovchain::markovchainFit(method = "map")` supports Bayesian
  estimation but only with uniform priors.
- \[markovchain package\]: Installed. `markovchainFit()` supports
  `method = "mle"`, `"laplace"` (add-1 smoothing), `"map"` (Bayesian).
  MAP returns `$estimate`, `$expectedValue`, `$standardError`,
  `$confidenceInterval`, `$logLikelihood`. Only supports uniform priors
  — no informative prior input.
- \[existing comparison functions\]: All TNA-specific.
  [`compare_networks()`](https://pak.dynasite.org/Saqrlab/reference/compare_networks.md)
  compares two TNA models (correlation, RMSE, edge weights).
  [`compare_network_estimation()`](https://pak.dynasite.org/Saqrlab/reference/compare_network_estimation.md)
  compares TNA model types (tna/atna/ctna/ftna) via sampling stability.
  [`compare_edge_recovery()`](https://pak.dynasite.org/Saqrlab/reference/compare_edge_recovery.md)
  does binary edge TP/FP/FN.
  [`cross_validate_tna()`](https://pak.dynasite.org/Saqrlab/reference/cross_validate_tna.md)
  cross-validates TNA model types. None work with `psych_network`
  objects from `build_network()`.
- \[estimate_network\]: Universal estimation engine.
  `estimate_network(data, method, params, scaling, threshold, level, id_col)`
  → `saqr_network` S3 class. Supports 6 built-in methods: relative,
  frequency, co_occurrence (transition); cor, pcor, glasso
  (association). `params` list is the composability key — stores
  method-specific args for replay by bootstrap/grid search.
  `build_network()` now delegates to `estimate_network()` internally and
  converts output to `psych_network` for backward compat via
  `.saqr_to_psych_network()`.
- 
- \[tabulate byrow\]: When using
  `pair_idx = (from_int - 1) * n_states + to_int` with
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html), the resulting
  flat vector is in row-major order. R’s
  [`matrix()`](https://rdrr.io/r/base/matrix.html) fills column-major by
  default, so `byrow = TRUE` is required to get correct mat\[from, to\]
  indexing. Without it, from/to get silently swapped.
- 
- \[count_transitions NA handling\]: When counting transitions, NAs must
  NOT be stripped before creating pairs. Instead, create consecutive
  pairs first (from = actions\[-n\], to = actions\[-1\]), then filter
  out pairs where either side is NA. Stripping NAs first creates false
  transitions by bridging gaps (e.g. A, NA, B → would falsely create
  A→B).
- \[frequencies vectorized\]: `frequencies()` now delegates to
  `.count_transitions()` which uses vectorized base R (matrix slicing +
  tabulate) for wide format and data.table for long format. Both are
  significantly faster than the old per-row `lapply` approach.
- \[bootstrap_network\]: Universal bootstrap for any estimator in
  `estimate_network()`. Two paths: transition methods use fast
  pre-computed per-sequence counts (single `tabulate` + `colSums` per
  iteration); association methods call the full estimator per iteration.
  Returns `saqr_bootstrap` S3 class with CIs, p-values, mean/SD,
  significant edges, pruned network. ~2.8x faster than
  [`tna::bootstrap()`](http://sonsoles.me/tna/reference/bootstrap.md) on
  2000 sequences / 500 iterations due to `colSums` on 2D matrix
  (C-level) vs `apply` on 3D array (R-level).
- \[bootstrap precomputation\]: Per-sequence count matrices are built as
  n_seq x (n_states^2) 2D matrices using combined flat index
  `(row - 1) * nbins + pair_idx` with single
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html). Each bootstrap
  iteration only needs `colSums(trans_2d[sampled_rows, ])` + lightweight
  post-processing (row-normalize for relative).
- 
- \[bootstrap summary\]: Long-format data frame built via
  [`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)
  for speed. Directed networks filter `weight != 0 & from != to`;
  undirected filter `weight != 0 & from < to` (upper triangle only).
- \[bootstrap refactor\]: Dropped long-format precompute helpers
  (bootstrap now requires wide format, errors with message to use
  `convert_sequence_format()`). Merged transition+co-occurrence wide
  precompute into single `.precompute_per_sequence_wide()`. Inlined
  `.postprocess_transition()` into vapply loop — builds matrix without
  dimnames since it’s flattened immediately. Extracted
  `.select_state_cols()` shared helper into estimators.R. File went from
  787 → 598 lines (~24% reduction). All 566 tests pass.
- \[build_network consolidation\]: `build_network()` is now the single
  entry point for all 6 estimation methods (was split between
  `build_network()` for association and `estimate_network()` for all).
  Returns `netobject` class (was `psych_network`/`saqr_network`).
  `estimate_network()` is now a deprecated wrapper that calls
  `build_network()`. `method` parameter is required (no default). Shared
  helpers (`.resolve_method_alias()`, `.apply_scaling()`,
  `.extract_edges_from_matrix()`, `.decompose_multilevel()`) live in
  estimate_network.R. S3 methods (print, plot, predictability) live in
  build_network.R.
- 
- 
- \[testthat edition 3 deprecation\]:
  [`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) emits a
  warning, and testthat edition 3 (set in DESCRIPTION
  `Config/testthat/edition: 3`) treats warnings as test failures. All
  calls to deprecated functions in tests must be wrapped in
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html).
- \[plot predictability conditional\]: `plot.netobject` only adds
  predictability pies for association methods
  (`method %in% c("glasso", "pcor", "cor")`), not transition methods.
  This avoids errors since transition methods don’t produce precision
  matrices.
- \[netobject \$data field\]:
  \`netobject\$data`stores the **cleaned** data from the estimator, not the raw input. Association methods: cleaned numeric matrix (NAs, non-numeric, zero-variance dropped by`.prepare_association_input()`). Transition methods: original data.frame as-is. Matrix input (e.g. pre-computed correlation): NULL (no row-level data). The`\$data\`
  field is not added to \`netobject_ml\` directly — each sub-network
  stores its own. Downstream code (\`permutation_test\`,
  \`bootstrap_network\`) can use \`\$data\` directly without
  re-cleaning.
- \[permutation_test\]: Universal permutation test for network
  comparison. Works with all 6 methods. Transition methods use fast
  precomputed counts (same as bootstrap). Association methods use
  optimized direct calls (see below). P-values:
  `(sum(|perm_diff| >= |obs_diff|) + 1) / (iter + 1)`. Effect size:
  `obs_diff / sd(perm_diffs)`. Running accumulators avoid storing full
  3D array. Supports paired mode (within-pair swaps) and p.adjust
  correction. S3 class: `saqr_permutation`.
- \[permutation_test optimization\]: Association path optimized for
  speed: (1) pre-clean pooled numeric matrix once before loop, (2)
  direct
  [`cor()`](https://rdrr.io/r/stats/cor.html)/[`solve()`](https://rdrr.io/r/base/solve.html)
  calls for cor/pcor (no estimator dispatch), (3) `glasso::glassopath()`
  with 50 lambda values + `.select_ebic_from_path()` for glasso.
  `glassopath()` is 14.5x faster than separate `glasso()` calls because
  it’s a single Fortran call for the full regularization path. Lambda
  path uses `.compute_lambda_path()` on pooled correlation matrix.
- \[permutation_test vs NCT\]: 100-dataset × 1000-iteration benchmark
  against `NetworkComparisonTest::NCT`: cor p-values **identical**
  (r=1.0, 100% agreement, 1.9x faster); glasso p-value r=0.9999, 99.9%
  agreement, 2.6x faster. Previous 98.6% cor disagreement was caused by
  RNG stream divergence (see sampling fix below).
- \[sample.int RNG divergence\]: `sample.int(n)` (full permutation)
  consumes `n` RNG values; `sample.int(n, k)` consumes only `k`. NCT
  uses partial draw `sample(1:n, k)`. Using `sample.int(n)` then taking
  first `k` causes RNG state to diverge from NCT after iteration 1,
  producing different permutations. Fix: use `sample.int(n_total, n_x)`
  to match NCT’s RNG consumption exactly. With same seed → identical
  results.
- \[permutation_test vs tna\]: 100-dataset simulation against
  [`tna::permutation_test`](http://sonsoles.me/tna/reference/permutation_test.md):
  identical p-values (r=1.0) and negligible effect size differences for
  both TNA (relative) and CNA (co_occurrence) methods. Confirms
  numerical correctness.
- \[glassopath\]:
  `glasso::glassopath(s, rholist, trace = 0, penalize.diagonal)` returns
  `$wi` as a 3D array (p × p × n_lambda). Use `.select_ebic_from_path()`
  to pick the best precision matrix via EBIC. Much faster than calling
  `glasso()` separately for each lambda because it warm-starts from the
  previous solution internally.
- \[fixed lambda rejected\]: Using a fixed regularization lambda (from
  the original network) across permutations drops significance agreement
  from 99% to 93%. Must re-select lambda via EBIC on each permutation
  for valid inference. The `glassopath()` approach preserves proper
  model selection while being fast.

### 2026-02-17

- 

### 2026-02-22

- \[GLLA factorial scaling\]: The GLLA Vandermonde basis `t^k` gives
  Taylor coefficients, not derivatives. Row k of the weight matrix must
  be multiplied by `k!` to get the actual k-th derivative. Without this,
  velocity (k=1) is correct (1!=1), but acceleration (k=2) is off by
  factor 2.
- \[.count_transitions_long time=NULL\]: Passing `time = NULL` to
  `.count_transitions_long()` causes an error because
  `NULL %in% names(dt)` returns `logical(0)` which can’t be used in
  `if()`. When calling from within `.matrices_from_raw_data()` (where
  data is already split by time), pass the actual `time_col` name
  instead of NULL.
- \[velocity_tna\]: Implemented as `velocity_tna()` in
  `R/velocity_tna.R`. Supports three methods: `"regression"` (default,
  OLS or beta), `"glla"`, `"finite_difference"`. Three input formats:
  list of matrices, `cograph::tna_windows()` output, raw data frame.
  Returns `tna_velocity` S3 class with print/summary/plot methods. Plot
  supports “network” (cograph::splot with green/red edges), “series”
  (matplot with regression fit lines), and “heatmap” (image with
  diverging palette). Regression returns edge_stats with slope, SE,
  t-value, p-value, R² per edge.
- \[as.vector column-major\]: `as.vector(matrix)` in R flattens
  column-major. When iterating over flattened edge indices j=1..n², the
  correct mapping is `row_idx = ((j-1) %% n) + 1`,
  `col_idx = ((j-1) %/% n) + 1`. Swapping `%%` and `%/%` transposes the
  result silently — a subtle bug that passes dimension checks but puts
  values in wrong cells.
- \[Ising model estimator\]: Implemented as `.estimator_ising()` in
  `R/estimators.R`. Uses nodewise L1-penalized logistic regression
  (glmnet) with EBIC model selection (IsingFit algorithm, van Borkulo et
  al. 2014). Parameters: `gamma` (EBIC, default 0.25), `nlambda`
  (default 100), `rule` (“AND”/“OR” symmetrization). Requires `glmnet`
  (in Suggests). Method name: `"ising"`, alias: `"isingfit"`. Returns
  undirected weighted adjacency + thresholds. Validated: exact numerical
  match with IsingFit::IsingFit across 50+ random configurations.
- \[Ising EBIC formula\]: IsingFit uses **nodewise** EBIC:
  `-2*loglik + k*log(n) + 2*gamma*k*log(p-1)` where p-1 is the number of
  predictors per regression. This differs from graphical EBIC (used in
  EBICglasso) which uses `4*gamma*k*log(p)`. Using the wrong formula
  selects different lambdas for ~20% of cases. IsingFit’s `min_sum`
  defaults to `-Inf` (no row filtering).
- \[Ising OR symmetrization\]: IsingFit’s OR rule uses simple
  `(W + t(W)) / 2` — averaging ALL values including zeros. NOT “average
  of nonzero values only”. The AND rule masks by `(W!=0) & (t(W)!=0)`
  first, then averages.
- 

### 2026-02-23

- \[temporal BFS convergence\]: Temporal BFS (earliest-arrival) can
  require multiple passes when edges are not in strict temporal DAG
  order (i.e., when earlier-onset edges can update arrivals that unlock
  later-onset edges already scanned). Using a `changed` flag with
  `while` loop bounded by `n_vertices` iterations guarantees
  convergence. For pure temporal DAGs, single pass suffices.
- \[backward temporal reachability\]: Implemented by reversing edges and
  negating times, then running forward BFS. The target vertex gets
  `start_time = -Inf`, which
  [`is.finite()`](https://rdrr.io/r/base/is.finite.html) returns FALSE
  for — so the target is already excluded from the “reachable” count. Do
  NOT subtract 1 again.
- \[sd() with single value\]: [`sd()`](https://rdrr.io/r/stats/sd.html)
  on a length-1 vector returns `NA`, not 0. When computing burstiness
  coefficient `B = (sigma - mu) / (sigma + mu)`, must handle
  `length(iet) < 2` by setting `sigma = 0` to avoid `NA` propagation.
- \[temporal density normalization\]:
  `density = sum(durations) / (n*(n-1) * time_span)` for directed
  networks. When a single edge scales proportionally with time_span
  (longer edge = longer span), density stays constant. To test density
  differences, vary the number of edges or add parallel edges within the
  same time span.
- \[temporal_network\]: Implemented in `R/temporal_network.R`. Input:
  edge list with (from, to, onset, terminus). Core algorithms: temporal
  BFS (earliest-arrival), time-varying degree, forward/backward
  reachability, temporal closeness, temporal betweenness,
  formation/dissolution, edge durations, IET + burstiness, temporal
  density, static snapshots. No new dependencies (igraph + ggplot2
  already in Imports). Returns `temporal_network` S3 class with
  print/summary/plot methods. 192 tests (including 20 tsna equivalence
  tests).
- 
- \[tsna tDegree is point-in-time\]: `tsna::tDegree()` computes degree
  at each time point instant (`onset <= t & terminus > t`), NOT interval
  overlap (`onset < bin_end & terminus > bin_start`). For
  integer-aligned times, both agree. For floating-point times, they
  differ — edges starting mid-bin are counted by interval overlap but
  not by the point-in-time snapshot. Use
  `onset <= bin_start & terminus > bin_start` to match tsna.
- \[tsna tEdgeFormation uses exact equality\]: `tsna::tEdgeFormation()`
  counts edges where `onset == t` (exact equality). This only matches
  our interval-based count (`onset >= bin_start & onset < bin_end`) when
  onset times are grid-aligned (integer). For floating-point times, tsna
  misses most edges since `onset == grid_point` is rarely true. Our
  interval approach is more general; only compare with tsna on
  integer-time datasets.
- \[tsna tEdgeDissolution convention\]: tsna records dissolution at the
  terminus time (the first moment the edge is gone). Uses
  `[bin_start, bin_end)` half-open interval, same as formation. Our
  original `(bin_start, bin_end]` convention was shifted by one bin.
- \[networkDynamic merges overlapping spells\]: When two edge spells for
  the same (tail, head) pair overlap in time, `networkDynamic` merges
  them into a single activity period. This changes total edge duration
  (union vs sum). Our implementation sums durations per pair without
  merging overlaps. Compare edge durations with tsna only on datasets
  without overlapping spells for the same pair.
- \[numeric vertex sorting\]: When vertex names are all numeric strings
  (“1”, “2”, …, “10”), [`sort()`](https://rdrr.io/r/base/sort.html)
  gives alphabetical order (“1”, “10”, “2”, …). Must detect numeric
  names and sort by
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) to match tsna’s
  integer vertex indexing.
- \[undirected temporal BFS\]: Temporal BFS only traverses from→to by
  default. For undirected networks, must add reverse edges (to→from with
  same onset/terminus) before running BFS, otherwise reachability is
  asymmetric. Both `.compute_all_temporal_paths()` and
  `temporal_paths()` need this fix.
- 
- 
- 
- \[igraph empty graph metrics\]: On empty graphs (0 edges):
  `reciprocity()` → NaN, `transitivity(type="global")` → NaN,
  `centr_clo()$centralization` → NaN, `mean_distance(unconnected=TRUE)`
  → NaN, `assortativity_degree()` → NaN, `constraint()` → NaN for all
  vertices. `page_rank()` → uniform 1/n. `hub_score()/hits_scores()$hub`
  → all 1. `coreness()` → all 0. `centr_degree()` → 0. `centr_betw()` →
  0.
- 
- \[sna vacuous transitivity\]:
  [`sna::gtrans()`](https://rdrr.io/pkg/sna/man/gtrans.html) returns 1
  when there are no two-paths (vacuously true).
  [`igraph::transitivity()`](https://r.igraph.org/reference/transitivity.html)
  returns NaN. For tsna parity, use 1 for the vacuous case.
- \[snapshot point-in-time parity\]: tsna’s `tSnaStats()` uses
  `network.collapse()` which snapshots edges active at a specific time
  point: `onset <= t & terminus > t`. Our original interval-overlap
  approach (`onset < bin_end & terminus > bin_start`) included edges not
  yet active at bin_start but starting within the bin. Changed to
  point-in-time +
  [`igraph::simplify()`](https://r.igraph.org/reference/simplify.html)
  for exact tsna match.
- 
- 
- \[centr_degree directed star\]: For a directed out-star
  (A→B,A→C,A→D,A→E), `centr_degree(g)` with default mode=“all” (total
  degree) gives centralization ~0.375 (not near 1). The out-degree
  centralization via `centr_degree(g, mode="out")` gives 0.8. Default
  mode averages in+out degrees which dilutes the star pattern.

### 2026-02-26

- \[simulate_data cholesky\]: When generating multivariate data from a
  correlation matrix, use `chol(R)` (upper Cholesky) and `Z %*% L` where
  Z is iid standard normal. The `.nearest_pd()` helper (eigenvalue
  clamping + unit diagonal normalization) ensures the matrix is always
  valid even after arbitrary off-diagonal assignment.
- \[simulate_data factor analysis n constraint\]:
  [`factanal()`](https://rdrr.io/r/stats/factanal.html) requires n \> p
  (observations \> variables). The generator enforces
  `n >= max(p + 10, 50)` to avoid convergence issues. With 2–4 factors ×
  3–4 items = 6–16 variables, minimum n is 50.
- \[simulate_data seed-driven variation\]: Using
  [`sample()`](https://rdrr.io/r/base/sample.html) and
  [`runif()`](https://rdrr.io/r/stats/Uniform.html) at the start of each
  generator to pick structural parameters (n, k, d, effect sizes) from
  ranges makes different seeds produce structurally different datasets.
  The seed also governs the data itself, so reproducibility is exact.

### 2026-02-27

- \[boot_glasso CS-coefficient\]: CS-coefficient matching bootnet
  requires three key details: (1) Each case-drop iteration randomly
  picks ONE drop proportion from the cs_drop levels (bootnet’s
  random-per-iteration approach, not fixed iterations per
  proportion). (2) Uses Pearson correlation with original centralities,
  NOT Spearman — bootnet’s `cor0()` is Pearson. (3) CS = max(drop_prop)
  where `mean(cors > 0.7) > 0.95` (proportion-above-threshold, not
  5th-percentile-based). With these three fixes, CS matches
  `bootnet::corStability()` within 0.15 on 100% of seeds tested (mean
  diff 0.006).
- \[boot_glasso input dispatch\]: `.prepare_association_input()` treats
  raw matrices (n x p) as correlation matrices (checks nrow == ncol).
  Must convert to data.frame first when passing raw observation data as
  a matrix.
- 
- \[boot_glasso edge diff scalability\]: Pairwise edge difference test
  is O(n_edges^2 \* iter). For p=30 (435 edges), that’s 94,830 pairs.
  The function returns NULL for n_edges \> 500 to avoid excessive
  computation. For large networks, use centrality difference tests
  instead.
- \[mlvar within-centering\]: Person-mean centering BOTH Y (outcome) and
  X (predictor) in panel VAR absorbs the random intercept, making OLS
  equivalent to the fixed-effects estimator (same as
  `lmer(..., (1|id))`). Must center both sides — centering only Y leaves
  the intercept in X.
- \[mlvar df correction\]: After within-centering, effective df =
  n_obs - d - n_subjects (not n_obs - d). The n_subjects penalty
  accounts for the absorbed person-specific intercepts. SE_correct =
  SE_ols \* sqrt(df_ols / df_correct).
- \[glasso symmetry\]: `glasso::glasso()` precision matrix (`$wi`) can
  have floating-point asymmetry ~1e-6. `.precision_to_pcor()` propagates
  this. Must force symmetry with `(pcor + t(pcor)) / 2` after
  conversion.
- \[mlvar VAR simulation\]: To simulate VAR(1) with known B:
  `y_t = mu + B %*% (y_{t-1} - mu) + noise`. Person-specific `mu`
  creates between-subject variation. Innovation noise through correlated
  Cholesky creates contemporaneous structure. Diagonal of B controls
  autoregression (0.2–0.5 range works well for ESM data).
- 
- \[mlvar vs mlVAR equivalence\]: Our `mlvar(standardize=FALSE)`
  produces **identical** temporal coefficients to
  `mlVAR::mlVAR(estimator="lmer", temporal="fixed", scale=FALSE)` — max
  diff ~1e-15 across 25 seeds. P-values differ by \<0.002 (our
  t-distribution vs mlVAR’s normal approximation). Residual correlations
  match r\>0.9999 (lmer REML vs OLS). Person mean correlations match
  r\>0.999. Significance agreement \>95% at alpha=0.05.
- \[mlVAR estimator=“lm”\]: The `"lm"` estimator requires
  `temporal="unique"` (no other option). It uses LSDV (person dummies)
  instead of centering — algebraically equivalent but gives different
  scaling of coefficients due to how between-person means are included.
  The `"lmer"` estimator with `temporal="fixed"` is closer to our
  centering approach.
- \[mlVAR contemporaneous\]: mlVAR uses unregularized nodewise
  regression on residuals + matrix inversion for partial correlations
  (no GLASSO). Our GLASSO-based approach produces sparser networks. The
  underlying residual correlation matrices match, but the final partial
  correlations differ by method design.
- \[mlVAR p-value method\]: mlVAR uses `2*(1-pnorm(|beta/se|))` (normal
  approximation). Our implementation uses
  `2*pt(|t|, df, lower.tail=FALSE)` (exact t-distribution). With 1000+
  observations, the difference is \<0.002.

### 2026-02-28

- \[proximity timeline\]: `temporal_network()` uses `time_interval`
  parameter (not `n_bins`) to control temporal resolution. The number of
  bins is derived as `ceiling((t_max - t_min) / time_interval)`.
  Pre-computed metric matrices (eigenvector, page_rank, etc.) are all
  n_vertices × n_bins with vertex_names as rownames.
- \[cmdscale on disconnected graphs\]:
  [`igraph::distances()`](https://r.igraph.org/reference/distances.html)
  returns Inf for unreachable vertex pairs. Must replace Inf with a
  large finite value (e.g. n_vertices) before passing to
  [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html),
  otherwise MDS fails. Wrap in
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html) for degenerate
  cases (e.g. zero-edge snapshots).
- \[ggplot2 legend auto-hide\]: Setting
  `theme(legend.position = "none")` hides the legend. For proximity
  timeline, auto-hide when n_vertices \> 15 and no vertex_group
  provided, since individual vertex colors become unreadable.

### 2026-03-05

- \[cograph::mcml naming conflict\]:
  [`cograph::mcml()`](https://sonsoles.me/cograph/reference/mcml.html)
  is a deprecated alias for `cluster_summary()`. Saqrlab’s MCML function
  must be named differently (`build_mcml`) to avoid masking when both
  packages are loaded.
  [`library(cograph)`](https://sonsoles.me/cograph/) masks Saqrlab’s
  export silently.
- \[build_network matrix input\]: `build_network()` with
  `method = "relative"` does NOT accept a pre-computed matrix — the
  `relative` estimator calls `.count_transitions()` which requires a
  data.frame. For pre-computed matrices, pass directly to
  [`cograph::cluster_summary()`](https://sonsoles.me/cograph/reference/cluster_summary.html)
  instead.
- \[cluster_summary type=“tna”\]:
  `cograph::cluster_summary(..., type = "tna")` row-normalizes the
  between-cluster matrix (diagonal = 0) and each within-cluster matrix.
  The `$between` and `$within` elements are tna-class objects with
  `$weights`, `$inits`, `$labels`, `$data`.
- \[cograph as_cograph edge list\]:
  `cograph::as_cograph(df, directed=TRUE)` accepts edge list data.frame
  with from/to/weight columns. `cograph::to_matrix(g)` converts to
  square matrix. Nodes in `g$nodes` have `label` and `name` columns.
  Custom columns (like `cluster`) can be added to `g$nodes` for
  downstream use.
- \[mcml recode vs cluster_summary\]:
  [`cograph::cluster_summary()`](https://sonsoles.me/cograph/reference/cluster_summary.html)
  aggregates a node-level matrix post-hoc (zeroes diagonal, averages
  blocks, renormalizes). The recode/filter approach estimates directly
  from sequences: recode states to cluster labels for between, filter
  non-member states to NA for within. These produce fundamentally
  different results — the recode approach preserves actual transition
  patterns including self-loops and correct normalization. Always use
  recode/filter when raw sequences are available.
- \[bootstrap_mcml recode\]: To bootstrap between-cluster network,
  recode all states to their cluster labels in the wide data, then call
  `bootstrap_network()` on the recoded data. For within-cluster, filter
  data so non-cluster states become NA (bootstrap_network handles NAs by
  treating them as sequence breaks).
- \[mcml TNA alignment\]: Between-cluster network includes self-loops
  (non-zero diagonal), matching how TNA counts all consecutive
  transitions. Do NOT collapse consecutive same-cluster states — user
  explicitly rejected this approach. The between matrix is a standard
  TNA model on recoded data.

### 2026-03-06

- \[gimme algorithm\]: GIMME (Group Iterative Multiple Model Estimation)
  uses uSEM: y(t) = A*y(t) + B*y(t-1) + e(t). A = contemporaneous
  effects (DAG), B = lagged effects. Path search via lavaan modification
  indices (MI). Group paths: MI significant for \>75% of converging
  subjects. Individual paths: Bonferroni-corrected MI per person.
  Pruning: remove group paths with non-significant z-values across \>50%
  of subjects.
- \[gimme lavaan syntax\]: Base syntax needs explicit variance (`~~`),
  intercept (`~1`), and exogenous covariance terms for all lagged
  variables. Without exogenous covariances, lavaan treats lagged vars as
  having zero correlation, breaking MI calculation. AR paths
  (`V1 ~ V1lag`) can be optional via `ar` parameter.
- \[gimme vs gimme package\]: **Exact equivalence** — path counts
  identical, standardized coefficients diff=0 across 100 randomly
  generated datasets (70 × 3-var, 30 × 4-var, random A/Phi/Psi, 8-15
  subjects, 80-150 obs). Three fixes were needed: (1) standardized betas
  via `lavInspect(fit, "std")$beta` (not
  `parameterEstimates()$est`), (2) standardized z-values via
  `standardizedSolution()$z` for pruning (not
  `parameterEstimates()$z`), (3) resume forward search after
  nonconv/stability removal, not just after z-pruning.
- 
- \[gimme MI extraction\]:
  [`lavaan::modindices()`](https://rdrr.io/pkg/lavaan/man/modificationIndices.html)
  returns data.frame with `lhs`, `op`, `rhs`, `mi` columns. Filter to
  `op == "~"` (regression paths only). Match candidate paths by
  `lhs_rhs` concatenation. MI \> 0 and model converges = eligible for
  selection.
- \[gimme testWeights\]: Stability check on standardized beta matrix.
  Extracts `lavInspect(fit, "std")$beta`, subsets to endo rows × coln
  columns (= c(lag_names, endo_names)). Tests eigenvalues of lagged
  block (cols 1:p) and contemporaneous block (cols p+1:2p). If any
  Re(eigenvalue) \>= 1 → unstable. In individual search, unstable models
  trigger path removal (pop last added) until stable, then z-prune, then
  optionally resume forward search excluding removed paths. Without
  this, coefficients can differ by up to ~0.05 from gimme.
- \[gimme standardized output\]: gimme stores standardized betas
  (`lavInspect(fit, "std")$beta`) in `path_est_mats`, rounded to 4
  digits. Using unstandardized `parameterEstimates()$est` produces
  different values. Must use standardized for exact match.
- \[gimme standardized z-values\]: gimme uses `standardizedSolution()$z`
  for z-pruning in individual search (via `return.zs`), NOT
  `parameterEstimates()$z`. The standardized z can differ substantially
  from unstandardized (e.g. 2.90 vs 2.61), flipping prune decisions near
  the Bonferroni cutoff.
- 

### 2026-03-10

- \[BuildHON tuple keys in R\]: Python uses tuples as dict keys. R has
  no hashable vector type. Solution: encode character vectors with
  `paste(x, collapse = "\x01")` (non-printable separator). Decode with
  `strsplit(key, "\x01", fixed = TRUE)[[1]]`. Use R environments (hash =
  TRUE) for O(1) lookup.
- 
- \[BuildHON KLD division by zero\]: In theory, `KLD(a, b)` could hit
  `p_a/0` when a target exists in extended source but not in base
  source. In practice, this can’t happen with the original algorithm
  because every occurrence of extended source (A,B)→C is also counted as
  base source (B)→C, so if extended passes min_freq, base also passes.
  Safe to return Inf for the theoretical case.
- \[BuildHON SequenceToNode\]: The pipe notation reverses the history:
  `("X","A","B")` → `"B|A.X"` (current state first, then history
  most-recent-first, dot-separated). Python pops from the end of the
  remaining seq; R equivalent is `rev(seq[-length(seq)])`.
- 
- \[pyHON vs pathpy3\]: pyHON builds variable-order HON (the novel
  contribution). pathpy3 builds fixed-order HON (you specify k).
  Different algorithms, different outputs. Validate against pyHON, not
  pathpy3.
- \[reticulate RETICULATE_PYTHON\]: On macOS with homebrew Python, must
  set `Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/bin/python3")`
  before calling `py_available()`. Without this, reticulate may fail to
  find Python even when it’s installed.
- \[BuildHON min_freq threshold bug\]: Python’s `BuildDistributions`
  zeros out `Count[Source][Target]` in place when below `min_freq`. Then
  `KLDThreshold` uses `sum(Count[ext].values())` which includes the
  zeros — so the total reflects only surviving counts. R must also zero
  out `count` in place (not just filter when computing distributions),
  otherwise `.hon_kld_threshold` uses unfiltered totals, producing lower
  thresholds and spurious higher-order rules. Found via group_regulation
  dataset where R had 86 edges vs Python’s 78.
- \[reticulate dollar assignment in testthat::test_file\]:
  `reticulate::py$x <- val` inside a helper called from
  [`testthat::test_file()`](https://testthat.r-lib.org/reference/test_file.html)
  causes “object ‘reticulate’ not found” error because the `$<-`
  operator on a namespace-qualified object requires the package to be
  attached (not just namespace-loaded). Fix: use
  `reticulate::py_set_attr(reticulate::py, "x", val)` instead.
- \[reticulate py_to_r in test_file\]: When
  [`testthat::test_file()`](https://testthat.r-lib.org/reference/test_file.html)
  runs tests, iterating a Python list object (class
  `python.builtin.list`) with
  [`lapply()`](https://rdrr.io/r/base/lapply.html) fails because
  auto-conversion is not triggered in the test environment. Fix: call
  `reticulate::py_to_r()` on the Python list first to get a proper R
  list, then iterate with
  [`lapply()`](https://rdrr.io/r/base/lapply.html).
- \[reticulate Python return types for R\]: Python tuples returned by
  reticulate are converted to R environments in some contexts (not
  indexable by \[\[1\]\]). Returning Python dicts instead (as
  `{'key': val}`) gives named R lists that are reliably accessible via
  `e[["key"]]`.

### 2026-03-11

- \[HON node naming\]: Changed from pipe notation (`"B|A.X"`) to
  readable arrow notation (`"X -> A -> B"`). The pipe notation placed
  current state first with reversed history (dot-separated). Arrow
  notation reads naturally left-to-right as the state sequence.
  `.hon_sequence_to_node()` is now simply
  `paste(seq, collapse = " -> ")`.
- \[HON edge columns unified\]: HON edges now match MOGen/HYPA pattern:
  `path` (full sequence), `from` (context/from-node), `to` (next state),
  `count`, `probability`. Previously had `from`/`to` (node names in pipe
  notation) and `weight` (probability only, no counts).
- \[HON count recovery\]: Raw counts for HON edges are recovered from
  the `count` environment by looking up
  `count_env[[source_key]][[next_state]]`, where `next_state` is always
  the last element of the decoded target tuple (invariant under
  rewiring).
- 

### 2026-03-16

- \[R partial matching of n_batch\]: `n_batch` as a function parameter
  is partially matched by `n=100` passed via `...`, because “n” is a
  prefix of “n_batch”. Fix: place `n_batch` AFTER `...` in the function
  signature — R only allows exact-name matching for parameters after
  `...`.
- \[simulate_data complexity\]: Generation-time edge cases (tiny_n,
  heavy_tailed, heteroscedastic, extreme_imbalance, multicollinear) must
  be checked inside each generator via `opts$complexity`. Post-injection
  cases (na, outliers, ties, duplicates, constant_col, all_na_col) are
  applied after generation in `.inject_complexity_post()` and are
  universal across all types.
- 
- \[%\|\|% in simulate_data.R\]: After removing utils.R (moved to
  Nestimate), %\|\|% must be defined locally in simulate_data.R. Placed
  at the bottom of the file — R closures look up `%||%` in the package
  namespace at call time (not definition time), so order within a file
  does not matter.
- 
- \[R matrix() fills column-by-column\]:
  `matrix(c(0, 0, 10, 10), nrow=2, ncol=2)` puts 0 in col 1 and 10 in
  col 2. Confusion: `matrix(c(0, 10, 0, 10), nrow=2, ncol=2)` gives BOTH
  columns \[0,10\] — both profiles get identical means. Always verify
  matrix layout with [`print()`](https://rdrr.io/r/base/print.html)
  before writing test assertions.
- 
- 
- \[all.equal relative vs absolute tolerance\]:
  `all.equal(x, y, tolerance=0.1)` uses RELATIVE tolerance (mean
  absolute diff / mean absolute target), not absolute. A max absolute
  diff of 0.04 can still fail this test if the target values are small.
  Use `max(abs(x - y)) < threshold` for absolute element-wise
  comparison.
- \[simulate_fa .nearest_pd not suitable for Sigma\]: `.nearest_pd()` in
  simulate_data.R normalises diagonals to 1 (produces correlation
  matrix). It is only suitable for `phi` (which is a correlation
  matrix), never for `Sigma` (covariance matrix). Use raw eigenvalue
  clamp
  `eig$vectors %*% diag(pmax(eig$values, 1e-6)) %*% t(eig$vectors)` for
  covariance matrices.
- \[simulate_fa Cholesky data generation\]: `L <- chol(Sigma)` gives
  upper Cholesky (t(L) %*% L = Sigma). Generate with `Z %*% L` where Z
  is n×p iid N(0,1). Then Cov(X) = t(L) %*% L = Sigma. Correct and fast.
- \[list(data, params) return pattern\]: Returning list(data, params)
  instead of data.frame + attr() makes ground truth accessible to any
  consumer (JSON, CLI, Python, JS) without R-specific knowledge. attr()
  is R-internal and invisible to jsonlite::toJSON().
