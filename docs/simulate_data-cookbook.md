# `simulate_data()` Cookbook

A task-oriented reference for AI assistants generating test data with `simulate_data()`.

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
