# `simulate_data()` Cookbook

A task-oriented reference for AI assistants generating test data with `simulate_data()`.

---

## Standalone Usage

The file is fully self-contained. No package installation required — just source it:

```r
source("R/simulate_data.R")   # only dependency: base R + stats (always available)
```

All helpers (`%||%`, `.nearest_pd()`) and generators are defined in the file itself.

---

## Quick Start

```r
# Single clean dataset (reproducible, no edge cases)
d <- simulate_data("ttest", seed = 42)

# Single dataset with random edge cases injected
d <- simulate_data("ttest", seed = 42, complexity = "auto")

# All 7 types × 1000 datasets each (full test battery)
batches <- simulate_data("batch", seed = 1)
```

---

## Dataset Types

| `type` | Columns | Intended analysis | Extra `attr()` |
|---|---|---|---|
| `"ttest"` | `group` (factor), `score` (numeric) | `t.test(score ~ group, data = d)` | — |
| `"anova"` | `group` (factor), `score` (numeric) | `aov(score ~ group, data = d)` | — |
| `"correlation"` | `x1`…`xp` (numeric, p = 4–7) | `cor(d)`, `pairs(d)` | — |
| `"clusters"` | `x1`…`xd` (numeric), `true_cluster` (int) | `kmeans(d[,-ncol(d)], centers = max(d$true_cluster))` | — |
| `"factor_analysis"` | `x1`…`xp` (numeric, p = 6–16) | `factanal(d, factors = attr(d, "n_factors"))` | `n_factors`, `loadings` |
| `"prediction"` | `y`, `x1`–`x4` (numeric), `cat1`, `cat2` (factor) | `lm(y ~ ., data = d)` | — |
| `"mlvar"` | `id`, `day`, `beep` (int), `V1`…`Vd` (numeric) | `mlvar(d, vars = attr(d,"vars"), id="id", day="day", beep="beep")` | `true_temporal`, `true_contemporaneous`, `vars` |

All datasets also carry `attr(d, "type")` and `attr(d, "complexity")`.

---

## Testing Scenarios

Each scenario shows: the goal, the call, and what to assert.

---

### Test NA handling

**Goal:** Verify a function uses `na.rm`, imputes, or errors informatively on missing data.

```r
d <- simulate_data("ttest", seed = 1, complexity = "na")
stopifnot(anyNA(d$score))                # NAs were injected
t.test(score ~ group, data = d)         # should handle NAs or error clearly

# For multi-column types
d <- simulate_data("correlation", seed = 1, n = 100, complexity = "na")
stopifnot(anyNA(d))
cor(d, use = "pairwise.complete.obs")   # test na handling strategy
```

**What was injected:** 5–15% NAs spread across 50–100% of numeric columns.

---

### Test rank-based / tied-value handling

**Goal:** Verify functions that depend on ranks (Wilcoxon, Kruskal-Wallis, Spearman) handle ties correctly.

```r
d <- simulate_data("ttest", seed = 1, n = 100, complexity = "ties")
n_unique <- length(unique(d$score))
stopifnot(n_unique < 0.5 * nrow(d))    # many tied values confirmed

wilcox.test(score ~ group, data = d)   # check for ties warning/correction
```

**What was injected:** Numeric columns discretised to ~30% unique values.

---

### Test outlier robustness

**Goal:** Verify a function degrades gracefully or flags outliers when extreme values are present.

```r
d <- simulate_data("ttest", seed = 1, n = 200, complexity = "outliers")

# Robust check: values more than 5 MADs from median
med <- median(d$score, na.rm = TRUE)
mad_val <- mad(d$score, na.rm = TRUE)
stopifnot(any(abs(d$score - med) > 5 * mad_val, na.rm = TRUE))

# Compare robust vs non-robust estimator
mean(d$score)
median(d$score)
```

**What was injected:** 1–5 values set to 3–8 SDs from the column mean.

---

### Test small-sample behaviour

**Goal:** Verify error messages, warnings, and fallbacks when n is very small.

```r
d <- simulate_data("ttest", seed = 1, complexity = "tiny_n")
stopifnot(nrow(d) < 30)

# Expect warnings about low power, small df, etc.
tryCatch(t.test(score ~ group, data = d), warning = function(w) message("Warning: ", w))

# For factor analysis (n barely exceeds p)
d <- simulate_data("factor_analysis", seed = 1, complexity = "tiny_n")
stopifnot(nrow(d) < ncol(d) + 15)
```

**What was injected:** n forced to 6–20 (ttest/anova), p + 3–8 (factor_analysis), 10–25 (prediction), 3–8 subjects (mlvar).

---

### Test normality-assumption violations

**Goal:** Verify robustness tests, transformation checks, or non-parametric fallbacks.

```r
d <- simulate_data("ttest", seed = 1, n = 200, complexity = "heavy_tailed")

# Check kurtosis is elevated (t-distribution has heavier tails than normal)
library(moments)
stopifnot(kurtosis(d$score) > 3.5)

# Compare parametric vs non-parametric
t.test(score ~ group, data = d)
wilcox.test(score ~ group, data = d)
```

**What was injected:** Data generated from t(2–5) instead of Normal.

---

### Test heteroscedasticity handling

**Goal:** Verify whether a function applies Welch correction, Levene's test, or weighted approaches.

```r
d <- simulate_data("ttest", seed = 1, n = 150, complexity = "heteroscedastic")

# Verify unequal variances
vars <- tapply(d$score, d$group, var)
stopifnot(max(vars) / min(vars) > 4)

# Welch t-test handles it; pooled does not
t.test(score ~ group, data = d, var.equal = FALSE)  # correct
t.test(score ~ group, data = d, var.equal = TRUE)   # wrong assumption

# anova equivalent
car::leveneTest(score ~ group, data = d)
```

**What was injected:** One group's SD scaled up by 3–8×.

---

### Test extreme group imbalance

**Goal:** Verify that functions handle near-empty minority groups without silent errors.

```r
d <- simulate_data("ttest", seed = 1, n = 100, complexity = "extreme_imbalance")

tab <- table(d$group)
stopifnot(min(tab) < 0.15 * nrow(d))   # minority group < 15% of n

# Many tests break silently or warn with tiny groups
t.test(score ~ group, data = d)

# Check both groups actually exist
stopifnot(length(tab) == 2L)
```

**What was injected:** Balance ratio of 3–10% / 90–97%.

---

### Test zero-variance column handling

**Goal:** Verify that functions detect or skip constant columns without division-by-zero errors.

```r
# One constant column
d <- simulate_data("correlation", seed = 1, n = 50, n_vars = 5,
                   complexity = "constant_col")
col_vars <- sapply(d, var, na.rm = TRUE)
stopifnot(any(col_vars == 0 | is.na(col_vars)))

cor(d)                    # should produce NaN column or error
cor(d[, col_vars > 0])   # safe: drop constant columns first

# All-NA column
d <- simulate_data("correlation", seed = 1, n = 50, n_vars = 5,
                   complexity = "all_na_col")
stopifnot(any(sapply(d, function(v) all(is.na(v)))))
```

**What was injected:** One secondary numeric column set to its mean (constant) or to all-NA.

---

### Test multicollinearity handling

**Goal:** Verify that regression/factor analysis/network methods handle near-singular matrices.

```r
d <- simulate_data("correlation", seed = 1, n = 100, complexity = "multicollinear")

R <- cor(d)
diag(R) <- 0
stopifnot(max(abs(R)) > 0.90)    # near-perfect correlation confirmed

# Regularised methods should handle this; unregularised may fail
solve(cor(d))                         # may be near-singular
glasso::glasso(cor(d), rho = 0.1)    # regularisation handles it

# For regression
d <- simulate_data("prediction", seed = 1, n = 100, complexity = "multicollinear")
fit <- lm(y ~ ., data = d)
car::vif(fit)                         # expect high VIF for x1/x2
```

**What was injected:** One column replaced with another + 5% noise (r ≈ 0.99).

---

### Test duplicate row handling

**Goal:** Verify that functions do not silently inflate n or produce wrong degrees of freedom.

```r
d <- simulate_data("ttest", seed = 1, n = 80, complexity = "duplicates")
stopifnot(nrow(d) > 80)              # rows were appended

n_unique_rows <- nrow(unique(d))
stopifnot(n_unique_rows < nrow(d))   # duplicates confirmed

# Check whether function accounts for effective n
t.test(score ~ group, data = d)      # inflated n → overconfident p-value
t.test(score ~ group, data = unique(d))   # deduplicated baseline
```

**What was injected:** 10–25% of rows duplicated and appended.

---

### Combine multiple violations

**Goal:** Stress-test a function against realistic messy real-world data.

```r
# Explicit combination
d <- simulate_data("ttest", seed = 1, n = 150,
                   complexity = c("na", "outliers", "heavy_tailed"))

# Random mix (driven by seed — fully reproducible)
d <- simulate_data("ttest", seed = 42, complexity = "auto")
cat("Applied cases:", attr(d, "complexity"), "\n")

# Typical real-world stress test for a regression function
d <- simulate_data("prediction", seed = 5,
                   complexity = c("na", "outliers", "multicollinear", "duplicates"))
```

**Note:** Generation-time cases (`tiny_n`, `heavy_tailed`, `heteroscedastic`, `extreme_imbalance`, `multicollinear`) affect data shape. Post-injection cases (`na`, `outliers`, `ties`, `duplicates`, `constant_col`, `all_na_col`) are applied after and can be combined freely.

---

## Batch Generation

### One type, many datasets

```r
# 100 ttest datasets with auto complexity
batch <- simulate_data("ttest", seed = 1, n_batch = 100)
length(batch)          # 100
class(batch)           # "sim_batch_type" "list"

# Access individual datasets
d1 <- batch[[1]]
attr(d1, "seed")       # exact seed used (integer)
attr(d1, "batch_id")   # 1
attr(d1, "complexity") # character vector of applied cases

# Run a test across the batch
p_vals <- vapply(batch, function(d) {
  tryCatch(
    t.test(score ~ group, data = d)$p.value,
    error = function(e) NA_real_
  )
}, numeric(1L))

mean(p_vals < 0.05, na.rm = TRUE)   # proportion significant
mean(is.na(p_vals))                 # proportion erroring out
```

### All types at once

```r
# 7 types × 1000 datasets (default)
all_batches <- simulate_data("batch", seed = 1)
names(all_batches)     # "ttest" "anova" "correlation" ...
length(all_batches$ttest)   # 1000

# Custom batch size
all_batches <- simulate_data("batch", seed = 1, n_batch = 200)

# Iterate over all types
results <- lapply(all_batches, function(type_batch) {
  vapply(type_batch, function(d) nrow(d), integer(1L))
})
```

### Reproducibility guarantee

```r
# Same seed always gives the same datasets
b1 <- simulate_data("ttest", seed = 42, n_batch = 5)
b2 <- simulate_data("ttest", seed = 42, n_batch = 5)
identical(b1, b2)   # TRUE

# Dataset i from a single-type batch matches dataset i in a full batch
single  <- simulate_data("ttest", seed = 42, n_batch = 10)
full    <- simulate_data("batch", seed = 42, n_batch = 10)
identical(single[[1]][, c("group","score")],
          full$ttest[[1]][, c("group","score")])   # TRUE
```

---

## What to Check

Quick assertion patterns by scenario:

| Scenario | Assert |
|---|---|
| NA injection | `anyNA(d)` or `anyNA(d$score)` |
| Ties | `length(unique(d$score)) < 0.5 * nrow(d)` |
| Outliers | `any(abs(d$score - median(d$score)) > 5 * mad(d$score), na.rm=TRUE)` |
| Tiny n | `nrow(d) < 30` |
| Heavy tails | `kurtosis(d$score) > 3.5` (requires `moments`) |
| Heteroscedastic | `max(tapply(d$score, d$group, var)) / min(...) > 4` |
| Extreme imbalance | `min(table(d$group)) < 0.15 * nrow(d)` |
| Constant column | `any(sapply(d, var, na.rm=TRUE) == 0, na.rm=TRUE)` |
| All-NA column | `any(sapply(d, function(v) all(is.na(v))))` |
| Multicollinear | `max(abs(cor(d)[upper.tri(cor(d))])) > 0.90` |
| Duplicates | `nrow(d) > nrow(unique(d))` |
| Applied cases | `attr(d, "complexity")` returns character vector |

---

## Parameter Reference

| Parameter | Type | Values | Default | Notes |
|---|---|---|---|---|
| `type` | character | `"ttest"`, `"anova"`, `"correlation"`, `"clusters"`, `"factor_analysis"`, `"prediction"`, `"mlvar"`, `"batch"` | required | `"batch"` generates all 7 types |
| `seed` | integer or NULL | any integer | `NULL` | Controls both structure and data |
| `complexity` | character | `"clean"`, `"auto"`, or vector of case names | `"clean"` | `"auto"` defaults to `"clean"` for non-batch single datasets |
| `...` | named args | `n`, `n_groups`, `n_vars`, `n_clusters`, `n_dims`, `n_factors`, `items_per_factor`, `n_cat1`, `n_cat2`, `effect_size`, `n_subjects`, `d`, `n_days`, `beeps_per_day` | — | Structural overrides; passed to generator |
| `n_batch` | integer or NULL | any positive integer | `NULL` (1000 for `"batch"`) | **Must be named explicitly** — do not use `n=` |

### Complexity cases

| Case | Kind | Types affected | What it does |
|---|---|---|---|
| `"na"` | post | all | 5–15% NAs in 50–100% of numeric columns |
| `"outliers"` | post | all | 1–5 values set to 3–8 SDs from column mean |
| `"ties"` | post | all | Numeric cols discretised to ~30% unique values |
| `"duplicates"` | post | all | 10–25% of rows duplicated and appended |
| `"constant_col"` | post | all (≥2 numeric cols) | One secondary numeric col set to constant |
| `"all_na_col"` | post | all (≥2 numeric cols) | One secondary numeric col set to all-NA |
| `"tiny_n"` | gen | all | Very small n (type-specific minimum) |
| `"heavy_tailed"` | gen | all | t(2–5) distribution instead of Normal |
| `"heteroscedastic"` | gen | ttest, anova | Group SDs differ by 3–8× |
| `"extreme_imbalance"` | gen | ttest, anova | 3–10% / 90–97% group size split |
| `"multicollinear"` | gen | correlation, prediction, factor_analysis | Near-duplicate column added (r ≈ 0.99) |

**gen** = generation-time (affects data shape/distribution). **post** = post-injection (applied after generation, combinable freely).

---

## Common Patterns for Automated Testing

```r
# Pattern 1: Parametric sweep over 100 seeds, clean data
p_vals <- vapply(1:100, function(s) {
  d <- simulate_data("ttest", seed = s)
  t.test(score ~ group, data = d)$p.value
}, numeric(1L))

# Pattern 2: Full edge-case battery for a single function
cases <- c("na", "outliers", "ties", "tiny_n", "heavy_tailed",
           "heteroscedastic", "extreme_imbalance", "constant_col",
           "duplicates", "multicollinear")

results <- lapply(cases, function(case) {
  d <- simulate_data("ttest", seed = 1, n = 100, complexity = case)
  tryCatch(
    list(ok = TRUE,  result = t.test(score ~ group, data = d)),
    error   = function(e) list(ok = FALSE, error = conditionMessage(e)),
    warning = function(w) list(ok = TRUE,  warning = conditionMessage(w))
  )
})
names(results) <- cases

# Pattern 3: Pre-generated reproducible test suite
test_suite <- simulate_data("batch", seed = 42, n_batch = 50)
# Save once, reuse across test runs — attr(d, "seed") allows exact replay
saveRDS(test_suite, "tests/fixtures/test_suite.rds")
```

---

## Parameter-Recovery Functions

Four companion functions generate data with **fully specified ground-truth parameters** so that statistical models can be tested for correct recovery.

```r
source("R/simulate_latent.R")   # only dependency: base R + stats
```

All four functions return `list(data, params)` — both the observations and the ground truth are directly accessible fields, with no R-specific `attr()` needed.

```r
r <- simulate_fa(loadings, n = 500, seed = 1)
r$data    # data.frame of observations
r$params  # ground-truth parameters — plain R list, JSON-serializable

# Export for JavaScript / CLI consumption
jsonlite::write_json(r, "output.json", auto_unbox = TRUE)
```

---

### `simulate_lpa()` — Latent Profile Analysis

Generates continuous-indicator data where each observation belongs to one of K latent profiles.

```r
simulate_lpa(means, sds, props, n, seed = NULL)
```

| Argument | Type | Description |
|---|---|---|
| `means` | matrix (n_vars × n_profiles) | Column means per profile per variable |
| `sds` | matrix OR scalar | SDs per variable per profile (recycled if scalar) |
| `props` | numeric vector length K | Mixing proportions (normalised internally) |
| `n` | integer | Total sample size |
| `seed` | integer or NULL | Reproducibility |

**Returns:** `list(data, params)`.
- `$data`: `data.frame` with columns `y1`…`yp` + `true_profile` (integer).
- `$params`: `list(means, sds, props)`.

```r
means <- matrix(c(0, 0, 10, 10), nrow = 2, ncol = 2)  # profile 1: (0,0), profile 2: (10,10)
r <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 500, seed = 1)

table(r$data$true_profile)                          # should be ~250 / 250
tapply(r$data$y1, r$data$true_profile, mean)        # should be near 0 and 10

# Recovery check
grp_means <- vapply(1:2, function(k) {
  colMeans(r$data[r$data$true_profile == k, c("y1","y2")])
}, numeric(2))
stopifnot(all(abs(grp_means - r$params$means) < 0.5))
```

---

### `simulate_lca()` — Latent Class Analysis

Generates binary-indicator data where each observation belongs to one of K latent classes.

```r
simulate_lca(item_probs, class_probs, n, seed = NULL)
```

| Argument | Type | Description |
|---|---|---|
| `item_probs` | matrix (n_items × n_classes) | P(item = 1 \| class) for each item × class |
| `class_probs` | numeric vector length K | Class mixing proportions |
| `n` | integer | Total sample size |
| `seed` | integer or NULL | Reproducibility |

**Returns:** `list(data, params)`.
- `$data`: `data.frame` with binary columns `item1`…`itemm` + `true_class` (integer).
- `$params`: `list(item_probs, class_probs)`.

```r
item_probs <- matrix(c(0.9, 0.1, 0.9, 0.1,
                       0.1, 0.9, 0.1, 0.9), nrow = 4, ncol = 2)
r <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                  n = 1000, seed = 1)

table(r$data$true_class)                              # should be ~500 / 500
tapply(r$data$item1, r$data$true_class, mean)         # should be near 0.9 and 0.1

# Recovery check
cl1 <- r$data[r$data$true_class == 1L, ]
stopifnot(abs(mean(cl1$item1) - r$params$item_probs[1, 1]) < 0.05)
```

---

### `simulate_regression()` — Linear Regression

Generates a regression dataset where `lm(y ~ ., data = r$data)` recovers the true coefficients at large n.

```r
simulate_regression(coefs, predictor_sds, error_sd, n, seed = NULL)
```

| Argument | Type | Description |
|---|---|---|
| `coefs` | named numeric vector | True coefficients including `"(Intercept)"` |
| `predictor_sds` | named numeric vector | SD of each predictor (generated as N(0, sd)) |
| `error_sd` | positive numeric | Residual SD |
| `n` | integer | Sample size |
| `seed` | integer or NULL | Reproducibility |

**Returns:** `list(data, params)`.
- `$data`: `data.frame` with `y` + one column per predictor.
- `$params`: `list(coefs, predictor_sds, error_sd)`.

```r
coefs <- c("(Intercept)" = 2, x1 = 3, x2 = -1, x3 = 0.5)
r <- simulate_regression(
  coefs         = coefs,
  predictor_sds = c(x1 = 1, x2 = 1, x3 = 1),
  error_sd      = 0.5, n = 1000, seed = 42
)

fit <- lm(y ~ ., data = r$data)
est <- coef(fit)
stopifnot(all(abs(est - r$params$coefs[names(est)]) < 0.2))
```

---

### `simulate_fa()` — Factor Analysis

Generates multivariate normal data from an explicit factor model. The model-implied covariance matrix `Sigma = loadings %*% phi %*% t(loadings) + diag(psi)` is returned directly for cross-checking.

```r
simulate_fa(loadings, phi = NULL, psi = NULL, n, seed = NULL)
```

| Argument | Type | Description |
|---|---|---|
| `loadings` | matrix (p × m) | Factor loading matrix. Must be a matrix. Zero = no path. |
| `phi` | matrix (m × m) | Factor correlation matrix. Defaults to `diag(m)` (orthogonal). |
| `psi` | numeric vector (length p) | Unique variances. Defaults to `1 - diag(L %*% phi %*% t(L))`. |
| `n` | integer | Sample size |
| `seed` | integer or NULL | Reproducibility |

**Returns:** `list(data, params)`.
- `$data`: `data.frame` with columns `y1`…`yp`.
- `$params`: `list(loadings, phi, psi, sigma_implied)`.

`sigma_implied` is the exact `p × p` covariance matrix the data were drawn from.

```r
# Two orthogonal factors, 3 items each
# matrix() fills column-by-column (default byrow = FALSE):
# col 1 (F1): 0.8, 0.7, 0.6, 0, 0, 0 — col 2 (F2): 0, 0, 0, 0.8, 0.7, 0.6
loadings <- matrix(c(0.8, 0.7, 0.6, 0,   0,   0,
                     0,   0,   0,   0.8, 0.7, 0.6), nrow = 6, ncol = 2)

r <- simulate_fa(loadings = loadings, n = 500, seed = 1)
r$data              # y1..y6
r$params$loadings   # input loadings (exact)
r$params$phi        # diag(2) — orthogonal
r$params$psi        # unique variances
r$params$sigma_implied  # 6×6 implied covariance

# Recovery: sample cov ≈ sigma_implied at large n
r2 <- simulate_fa(loadings = loadings, n = 5000, seed = 1)
max(abs(cov(r2$data) - r2$params$sigma_implied))  # should be < 0.1

# Oblique model (correlated factors)
phi <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
r3  <- simulate_fa(loadings = loadings, phi = phi, n = 500, seed = 1)
r3$params$phi   # factor correlation matrix

# JS export from CLI
# Rscript -e "source('R/simulate_latent.R'); r <- simulate_fa(loadings, n=500, seed=1);
#             jsonlite::write_json(r, 'fa_data.json', auto_unbox=TRUE)"
```

---

### `simulate_seq_clusters()` — Sequence Cluster Analysis

Generates wide-format Markov chain sequences where each row belongs to one of K clusters, each driven by its own transition matrix. Designed to test sequence clustering methods for correct recovery of cluster assignments and transition dynamics.

```r
simulate_seq_clusters(trans_list = NULL, props = NULL, n = 300, seq_length = 20,
                      init_probs = NULL, n_clusters = 3L, n_states = 10L,
                      states = NULL, seed = NULL)
```

| Argument | Type | Description |
|---|---|---|
| `trans_list` | list of matrices or `NULL` | K square row-stochastic transition matrices with matching `rownames`/`colnames`. `NULL` = auto-generate. |
| `props` | numeric vector length K | Mixing proportions (normalised internally). Defaults to equal mixing. |
| `n` | integer | Total number of sequences. Default 300. |
| `seq_length` | integer | Time points per sequence (columns `T1`…`T{seq_length}`). Default 20. |
| `init_probs` | vector, list, or `NULL` | Initial state distribution. Shared vector, per-cluster list, or `NULL` (uniform). |
| `n_clusters` | integer | Number of clusters when `trans_list = NULL`. Default 3. |
| `n_states` | integer | Number of states when `trans_list = NULL`. Default 10. |
| `states` | character vector or `NULL` | State names when `trans_list = NULL`. Defaults to `"S1"`, `"S2"`, … |
| `seed` | integer or NULL | Reproducibility. |

**Returns:** `list(data, params)`.
- `$data`: `data.frame` with columns `T1`…`T{seq_length}` (character state labels) + `true_cluster` (integer 1…K).
- `$params`: `list(trans_list, props, init_probs)` — the K matrices, normalised proportions, and per-cluster initial distributions.

#### Explicit mode — supply your own transition matrices

```r
# Two clusters: cluster 1 stays in A, cluster 2 stays in B
m1 <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE,
             dimnames = list(c("A","B"), c("A","B")))
m2 <- matrix(c(0.2, 0.8, 0.1, 0.9), nrow = 2, byrow = TRUE,
             dimnames = list(c("A","B"), c("A","B")))

r <- simulate_seq_clusters(trans_list = list(m1, m2), n = 200, seed = 1)

# Sequences and ground truth
head(r$data[, c("T1","T2","T3","true_cluster")])
table(r$data$true_cluster)   # ~100 per cluster (equal mixing by default)

# Ground truth parameters
r$params$trans_list[[1]]     # m1 (cluster 1 matrix)
r$params$props               # c(0.5, 0.5)
r$params$init_probs[[1]]     # uniform: c(A=0.5, B=0.5)
```

#### Auto mode — random matrices, control cluster count and state space size

```r
# 3 clusters, 5 states (default n_states = 10)
r <- simulate_seq_clusters(n = 300, n_clusters = 3, n_states = 5, seed = 42)

# State names default to S1..S5
unique(unlist(r$data[, grep("^T", names(r$data))]))  # "S1" "S2" "S3" "S4" "S5"

# Verify row-stochastic
vapply(r$params$trans_list, function(m) max(abs(rowSums(m) - 1)), numeric(1))  # all near 0
```

#### Unequal mixing proportions

```r
r <- simulate_seq_clusters(trans_list = list(m1, m2, m3),
                           props = c(0.5, 0.3, 0.2),
                           n = 3000, seed = 1)
r$params$props                          # c(0.5, 0.3, 0.2) — stored normalised
table(r$data$true_cluster) / 3000       # observed counts ≈ props at large n
```

#### Custom initial state distribution

```r
# Shared across all clusters
r <- simulate_seq_clusters(trans_list = list(m1, m2),
                           init_probs = c(A = 0.8, B = 0.2),
                           n = 200, seed = 1)

# Per-cluster (list of vectors)
r <- simulate_seq_clusters(
  trans_list = list(m1, m2),
  init_probs = list(c(A = 0.9, B = 0.1), c(A = 0.1, B = 0.9)),
  n = 200, seed = 1
)
r$params$init_probs[[1]]   # c(A=0.9, B=0.1)
r$params$init_probs[[2]]   # c(A=0.1, B=0.9)
```

#### Comparing recovered clusters to ground truth

```r
r <- simulate_seq_clusters(n = 300, n_clusters = 3, n_states = 10, seed = 42)

seq_cols <- grep("^T", names(r$data), value = TRUE)
# ... run your clustering method on r$data[, seq_cols] → my_clusters ...

# Agreement table
table(recovered = my_clusters, true = r$data$true_cluster)

# Accuracy (after resolving label permutation)
mean(my_clusters == r$data$true_cluster)
```

#### Validation errors

```r
# Non-square matrix
bad <- list(matrix(1:6, nrow = 2))
simulate_seq_clusters(trans_list = bad, n = 10)      # Error: "square"

# Rows don't sum to 1
bad_mat <- matrix(c(0.5, 0.5, 0.5, 0.1, 0.9, 0.1, 0.3, 0.3, 0.3),
                  nrow = 3, byrow = TRUE,
                  dimnames = list(c("A","B","C"), c("A","B","C")))
simulate_seq_clusters(trans_list = list(bad_mat), n = 10)  # Error: "sum to 1"

# Mixed dimensions
simulate_seq_clusters(trans_list = list(m2x2, m3x3), n = 10)  # Error: "same"

# No rownames
bad <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2)
simulate_seq_clusters(trans_list = list(bad), n = 10)      # Error: "rownames"
```

---

**What to check:**

| Goal | Assert |
|---|---|
| LPA profile means recovered | `abs(colMeans(r$data[r$data$true_profile==k, -ncol(r$data)]) - r$params$means[,k]) < tol` |
| LCA item probs recovered | `abs(tapply(r$data$item1, r$data$true_class, mean) - r$params$item_probs[1,]) < tol` |
| Regression coefs recovered | `abs(coef(lm(y~., data=r$data)) - r$params$coefs[names(coef(...))]) < tol` |
| FA covariance recovered | `max(abs(cov(r$data) - r$params$sigma_implied)) < 0.1` (large n) |
| Seq cluster assignment | `table(recovered = my_clusters, true = r$data$true_cluster)` |
| Seq cluster proportions | `abs(table(r$data$true_cluster)/n - r$params$props) < 0.05` (large n) |
| Seq transition matrices | `r$params$trans_list[[k]]` — the exact matrix used for cluster k |
| Ground truth accessible | `r$params$coefs` / `r$params$means` / `r$params$item_probs` / `r$params$sigma_implied` / `r$params$trans_list` |
