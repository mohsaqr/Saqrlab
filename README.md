# Saqrlab <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/mohsaqr/Saqrlab/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mohsaqr/Saqrlab/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**A modern laboratory for data simulation.**

Saqrlab generates synthetic datasets with **known ground-truth parameters** so you can test, validate, and benchmark statistical methods — does your estimator actually recover the truth you simulated? It spans classic designs (t-test, ANOVA, regression, factor analysis), modern data-generating processes (IRT, survival, multilevel, growth curves, hidden Markov, missing-data mechanisms), and a full Temporal Network Analysis (TNA) toolkit, all behind one consistent interface.

- **87 functions across 15 categories** — every one documented with a runnable example ([browse the HTML reference](reference/index.html)).
- **One return type, `saqr_sim`** — every simulator gives you `$data` (a tidy data frame) and `$params` (the true generating parameters).
- **A unified `simulate()` dispatcher**, scenario presets, and an automatic **`validate_recovery()`** scorer that checks whether a fitted method recovered the truth.
- **1,173 tests, 0 failures**; ships a cross-package fixture-contract guard so generated data stays reproducible.

> Network *estimation* (bootstrap, permutation, GLASSO, GIMME, MCML, …) lives in the sibling package [**Nestimate**](https://cran.r-project.org/package=Nestimate). Saqrlab is simulation-first.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mohsaqr/Saqrlab")
```

Core dependencies install automatically. A few simulators have optional, skip-guarded recovery checks that use `lme4`, `mirt`, `survival`, `igraph`, `network`, or `sna` — install them only if you want those extras.

## Quick start

### One front door for every simulator

```r
library(Saqrlab)

list_simulators()                                  # the full catalogue, as a tidy table

sim <- simulate("ttest", n_a = 50, n_b = 50, mean_a = 0, mean_b = 0.6, seed = 1)
sim                                                # a saqr_sim: prints type, dims, params
head(sim$data)                                     # the data
sim$params                                         # the ground truth
```

### Simulate truth, then check recovery

```r
# Simulate a regression with known coefficients ...
sim <- simulate_regression(
  coefs         = c("(Intercept)" = 2, x1 = 3, x2 = -1),
  predictor_sds = c(x1 = 1, x2 = 1),
  error_sd      = 1, n = 1000, seed = 1
)

# ... fit a model ...
fit <- lm(y ~ x1 + x2, data = sim$data)

# ... and score how well it recovered the truth.
estimates <- setNames(coef(fit), paste0("coefs.", names(coef(fit))))
recovery  <- validate_recovery(sim, estimates = estimates)
recovery            # per-parameter true vs estimate, error, within-tolerance
summary(recovery)   # one-row scorecard: % within tolerance, mean error
```

### Modern data-generating processes

```r
simulate_irt(n_persons = 500, n_items = 20, model = "2PL", seed = 1)   # item responses
simulate_survival(n = 300, censoring_rate = 0.3, seed = 1)            # time-to-event
simulate_mlm(n_clusters = 30, cluster_size = 20, icc = 0.1, seed = 1) # students-in-classes
simulate_hmm(n_sequences = 50, seq_length = 30, n_states = 2, seed = 1)

# Inject realistic missingness (MCAR / MAR / MNAR) as a first-class step
inject_missingness(sim$data, mechanism = "MAR", prop = 0.15, predictor = "x1", seed = 1)
```

### Ready-made experimental designs

```r
list_scenarios()                                   # named design recipes
runs  <- run_scenario("power_ttest", seed = 1)     # runs every case, returns saqr_sim objects
tidy_simulation_results(runs)                      # one tidy data frame across all cases
```

### Temporal Network Analysis

```r
library(tna)

model <- simulate_tna_network(n_states = 6, seed = 42)   # a fitted `tna` object
plot(model)
centralities(model)
```

## How it is organised: two simulation tiers

Saqrlab deliberately keeps **two complementary interfaces** — knowing which one you're using tells you what you get back:

| Tier | Functions | Returns | Use it to |
|------|-----------|---------|-----------|
| **Explicit-parameter** | `simulate_ttest()`, `simulate_irt()`, `simulate_mlm()`, … (and `simulate(type, …)`) | a **`saqr_sim`** object with `$data` + `$params` | recover known truth — you pass the parameters in and check they come back out |
| **Random-parameter** | `simulate_data(type, seed = i)` | a **bare `data.frame`** (params in attributes) | stress / robustness testing — the seed invents a structurally unique dataset |

`saqr_sim` objects behave like their data frame for convenience (`head()`, `dim()`, `[`, `as.data.frame()`), while keeping `$params` for the ground truth.

## Function catalogue

Each category is a self-contained HTML page with every function's signature and a **runnable example showing real output** (open [`reference/index.html`](reference/index.html) for the clickable index):

| Category | Funcs | What's inside |
|----------|:----:|---------------|
| [Statistical simulators](reference/01-statistical.html) | 5 | `simulate_ttest`, `simulate_anova`, `simulate_correlation`, `simulate_clusters`, `simulate_prediction` |
| [Latent-variable models](reference/02-latent.html) | 5 | `simulate_lpa`, `simulate_lca`, `simulate_regression`, `simulate_fa`, `simulate_seq_clusters` |
| [Longitudinal & multilevel](reference/03-longitudinal-multilevel.html) | 3 | `simulate_longitudinal` (VAR/ESM), `simulate_mlm`, `simulate_growth` |
| [Item Response Theory](reference/04-irt.html) | 1 | `simulate_irt` (1PL/2PL/3PL/GRM) |
| [Survival & hidden Markov](reference/05-survival-hmm.html) | 2 | `simulate_survival`, `simulate_hmm` |
| [Missing-data mechanisms](reference/06-missing-data.html) | 1 | `inject_missingness` (MCAR/MAR/MNAR) |
| [Random-parameter generation](reference/07-random-parameter.html) | 1 | `simulate_data` (15 types + complexity injection + batch) |
| [TNA simulation](reference/08-tna-simulation.html) | 16 | `simulate_tna_network(s)`, `simulate_group_tna_networks`, `simulate_htna/mlna/mtna`, `generate_probabilities`, `sample_tna`, … |
| [Sequences](reference/09-sequences.html) | 2 | `simulate_sequences`, `simulate_sequences_advanced` |
| [Networks & graphs](reference/10-networks-graphs.html) | 5 | `simulate_igraph`, `simulate_network`, `simulate_edge_list`, `simulate_onehot_data`, `simulate_long_data` |
| [Reference data & utilities](reference/11-reference-data.html) | 10 | `LEARNING_STATES`, `GLOBAL_NAMES`, `get_learning_states`, `select_states`, `validate_sim_params`, … |
| [Comparison & model fitting](reference/12-comparison-fitting.html) | 10 | `fit_network_model`, `compare_networks`, `compare_centralities`, `compare_edge_recovery`, `cross_validate_tna`, … |
| [Visualization](reference/13-visualization.html) | 3 | `plot_network_estimation`, `plot_sampling_distribution`, `plot_tna_comparison` |
| [Batch, grid & sampling](reference/14-batch-grid-sampling.html) | 14 | `generate_param_grid`, `run_grid_simulation`, `run_bootstrap_simulation`, `summarize_grid_results`, … |
| [Laboratory infrastructure](reference/15-lab-infrastructure.html) | 9 | `simulate`, `list_simulators`, `validate_recovery`, `list_scenarios`, `run_scenario`, `tidy_simulation_results`, `export_simulation`, `saqr_sim` |

## Learning states

Saqrlab ships 180+ learning-action verbs across 8 categories, used to give simulated TNA states human-readable names.

| Category | Examples |
|----------|----------|
| metacognitive | Plan, Monitor, Evaluate, Reflect, Regulate |
| cognitive | Read, Study, Analyze, Summarize, Connect |
| behavioral | Practice, Annotate, Research, Review, Revise |
| social | Collaborate, Discuss, Seek_help, Explain, Share |
| motivational | Focus, Persist, Explore, Create, Commit |
| affective | Enjoy, Appreciate, Value, Curious, Cope |
| group_regulation | Adapt, Cohesion, Consensus, Coregulate, Plan |
| lms | View, Access, Download, Submit, Navigate |

```r
list_learning_categories()                          # all categories, as a table
get_learning_states("metacognitive", n = 5, seed = 1)
select_states(10, primary_categories = "metacognitive", seed = 1)
```

## Quality & reproducibility

- **1,173 tests, 0 failures** (`Rscript -e 'devtools::test()'`); package R code passes `R CMD check` cleanly.
- Every simulator is reproducible: the same `seed` gives identical output.
- A **cross-package fixture-contract guard** pins the exact output of `simulate_data()` for the seeds that generate downstream JSON fixtures, so reproducibility can't silently drift.
- Every example in the HTML reference and in roxygen is executed — the outputs you see are real.

## Documentation

- **[HTML function reference](reference/index.html)** — one page per category, every function with a live example.
- [`FEATURES.md`](FEATURES.md) — concise feature list. [`NEWS.md`](NEWS.md) — changelog.
- `?Saqrlab` and `?<function>` in R for the manual pages.

## Citation

```
Saqr, M. (2025). Saqrlab: A Modern Laboratory for Data Simulation.
R package version 0.4.0. https://github.com/mohsaqr/Saqrlab
```

## License

MIT License — see [LICENSE](LICENSE).

## Author

**Mohammed Saqr** — [saqr@saqr.me](mailto:saqr@saqr.me). Contributions welcome via [GitHub](https://github.com/mohsaqr/Saqrlab).
