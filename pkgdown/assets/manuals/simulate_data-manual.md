# Saqrlab Data Simulation Manual

> Comprehensive reference for generating synthetic statistical datasets with known ground-truth parameters. Use these functions to test statistical methods, validate implementations, and benchmark analysis pipelines.

## Installation

```r
# Install from GitHub
devtools::install_github("mohsaqr/Saqrlab")
library(Saqrlab)
```

---

## Quick Start

Every explicit-parameter simulation function returns a **`saqr_sim`** object:

```r
r <- simulate_ttest(n_a = 50, n_b = 50, mean_a = 100, mean_b = 105, seed = 1)

r$data        # data.frame with the simulated data
r$params      # list of ground-truth generating parameters
r$type        # "ttest"
r$seed        # 1

# Convenience accessors:
r[1:5, ]          # first 5 rows of $data
as.data.frame(r)  # extracts $data
head(r)           # head of $data
nrow(r)           # number of rows
print(r)          # compact summary
summary(r)        # full summary with data + params overview
```

The `$params` list always contains the exact generating parameters, so you can verify parameter recovery:

```r
r <- simulate_anova(n = 500, means = c(10, 15, 20), sds = 2, seed = 42)
# True values:
r$params$means         # c(G1 = 10, G2 = 15, G3 = 20)
r$params$eta_squared   # true population eta-squared

# Recovered values:
fit <- aov(score ~ group, data = r$data)
summary(fit)
```

---

## Function Reference

### 1. simulate_ttest() --- Two-Group Comparison

Generates data for an independent-samples t-test with fully specified population parameters.

**Signature:**

```r
simulate_ttest(
  n_a,                     # Integer. Sample size for group A. Min: 2.
  n_b,                     # Integer. Sample size for group B. Min: 2.
  mean_a,                  # Numeric. Population mean for group A.
  mean_b,                  # Numeric. Population mean for group B.
  sd_a     = 1,            # Numeric (> 0). SD for group A.
  sd_b     = sd_a,         # Numeric (> 0). SD for group B. Defaults to sd_a.
  labels   = c("A", "B"),  # Character[2]. Group labels.
  seed     = NULL           # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `group` | factor | Group membership (levels from `labels`) |
| `score` | numeric | Observed scores |

| `$params` fields | Type | Description |
|---|---|---|
| `mean_a` | numeric | Population mean for group A |
| `mean_b` | numeric | Population mean for group B |
| `sd_a` | numeric | Population SD for group A |
| `sd_b` | numeric | Population SD for group B |
| `n_a` | integer | Sample size for group A |
| `n_b` | integer | Sample size for group B |
| `cohens_d` | numeric | True Cohen's d = (mean_b - mean_a) / pooled_sd |

**Examples:**

```r
# Basic two-group comparison
r <- simulate_ttest(n_a = 30, n_b = 30, mean_a = 50, mean_b = 55, seed = 1)
t.test(score ~ group, data = r$data)
r$params$cohens_d  # true effect size

# Unequal variance (Welch's t-test scenario)
r <- simulate_ttest(n_a = 50, n_b = 50, mean_a = 0, mean_b = 1,
                     sd_a = 1, sd_b = 3, seed = 42)
t.test(score ~ group, data = r$data, var.equal = FALSE)

# Large sample for power analysis
r <- simulate_ttest(n_a = 200, n_b = 200, mean_a = 100, mean_b = 103,
                     sd_a = 10, seed = 7)

# Custom labels
r <- simulate_ttest(n_a = 40, n_b = 40, mean_a = 5, mean_b = 7,
                     labels = c("control", "treatment"), seed = 99)
```

---

### 2. simulate_anova() --- Multi-Group Comparison

Generates data for a one-way ANOVA with fully specified group means and standard deviations.

**Signature:**

```r
simulate_anova(
  n,                  # Integer or integer vector. Per-group size (scalar) or
                      #   per-group sizes (vector of length k). Min per group: 2.
  means,              # Numeric vector. Population means (length k, k >= 2).
  sds    = 1,         # Numeric scalar or vector[k]. Group SDs. Scalar recycled.
  labels = NULL,      # Character vector[k] or NULL. Group labels.
                      #   NULL -> "G1", "G2", ...
  seed   = NULL        # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `group` | factor | Group membership |
| `score` | numeric | Observed scores |

| `$params` fields | Type | Description |
|---|---|---|
| `means` | named numeric vector | Population means per group |
| `sds` | named numeric vector | Population SDs per group |
| `n` | named integer vector | Per-group sample sizes |
| `labels` | character vector | Group labels |
| `eta_squared` | numeric | True population eta-squared |

**Examples:**

```r
# Equal groups, equal variance
r <- simulate_anova(n = 50, means = c(10, 12, 15), seed = 1)
summary(aov(score ~ group, data = r$data))
r$params$eta_squared

# Unequal groups, heteroscedastic
r <- simulate_anova(n = c(100, 50, 30), means = c(5, 8, 12),
                     sds = c(1, 2, 4), seed = 42)

# Custom labels
r <- simulate_anova(n = 40, means = c(60, 70, 75, 80),
                     labels = c("placebo", "low", "medium", "high"), seed = 7)

# Two groups (equivalent to paired-down t-test)
r <- simulate_anova(n = 100, means = c(0, 0.5), sds = 1, seed = 3)
```

---

### 3. simulate_correlation() --- Correlated Multivariate Data

Generates multivariate normal data from an explicit correlation or covariance matrix.

**Signature:**

```r
simulate_correlation(
  n,                    # Integer. Sample size. Min: 3.
  sigma,                # Numeric matrix (p x p). Correlation or covariance matrix.
                        #   Must be symmetric. Near-PD corrected internally.
  means     = NULL,     # Numeric vector[p] or NULL. Means. NULL -> all zeros.
  var_names = NULL,     # Character vector[p] or NULL. NULL -> "x1", "x2", ...
  seed      = NULL       # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `x1`, `x2`, ..., `xp` | numeric | Simulated variables (or custom `var_names`) |

| `$params` fields | Type | Description |
|---|---|---|
| `sigma` | numeric matrix | The input correlation/covariance matrix |
| `means` | named numeric vector | Population means |
| `is_correlation` | logical | TRUE if diagonal of `sigma` is all 1s |

**Examples:**

```r
# Correlation matrix
R <- matrix(c(1.0, 0.7, 0.3,
              0.7, 1.0, 0.5,
              0.3, 0.5, 1.0), nrow = 3)
r <- simulate_correlation(n = 500, sigma = R, seed = 1)
cor(r$data)            # should approximate R
r$params$sigma         # exact truth

# With means and custom names
r <- simulate_correlation(
  n = 1000, sigma = R,
  means = c(100, 50, 25),
  var_names = c("IQ", "GPA", "income"),
  seed = 42
)
colMeans(r$data)       # should approximate c(100, 50, 25)

# Covariance matrix (non-unit diagonal)
S <- matrix(c(4, 1.5, 1.5, 9), nrow = 2)
r <- simulate_correlation(n = 200, sigma = S, seed = 7)
r$params$is_correlation  # FALSE
cov(r$data)              # should approximate S
```

---

### 4. simulate_clusters() --- Cluster / Mixture Data

Generates data from a mixture of Gaussians with known cluster centers, standard deviations, and true membership labels.

**Signature:**

```r
simulate_clusters(
  n,                 # Integer (total) or integer vector[k] (per-cluster sizes).
  centers,           # Numeric matrix (k x d). Cluster centroids.
  sds   = 1,         # Scalar, vector[k], or matrix (k x d). SDs. Scalar recycled.
  props = NULL,       # Numeric vector[k] or NULL. Mixing proportions.
                      #   NULL -> equal mixing. Only used when n is scalar.
  seed  = NULL         # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `x1`, `x2`, ..., `xd` | numeric | Observed dimensions |
| `true_cluster` | integer | True cluster membership (1 to k) |

| `$params` fields | Type | Description |
|---|---|---|
| `centers` | numeric matrix (k x d) | Cluster centroids |
| `sds` | numeric matrix (k x d) | Full SD matrix |
| `props` | numeric vector | Normalised mixing proportions |
| `n` | named integer vector | Per-cluster sample sizes |

**Examples:**

```r
# Three 2D clusters
centers <- matrix(c(0, 0,
                    5, 5,
                    10, 0), nrow = 3, byrow = TRUE)
r <- simulate_clusters(n = 300, centers = centers, seed = 1)
plot(x2 ~ x1, data = r$data, col = r$data$true_cluster)

# Unequal cluster sizes
r <- simulate_clusters(n = c(200, 50, 50), centers = centers,
                        sds = c(0.5, 1.0, 2.0), seed = 42)

# Imbalanced proportions (total n allocated by props)
r <- simulate_clusters(n = 1000, centers = centers,
                        props = c(0.7, 0.2, 0.1), seed = 7)
table(r$data$true_cluster)

# Per-cluster, per-dimension SDs
sd_mat <- matrix(c(0.5, 1.0,   # cluster 1: tight in x1, wide in x2
                    2.0, 0.5,   # cluster 2: wide in x1, tight in x2
                    1.0, 1.0),  # cluster 3: uniform
                  nrow = 3, byrow = TRUE)
r <- simulate_clusters(n = 300, centers = centers, sds = sd_mat, seed = 3)

# High-dimensional (5D, 4 clusters)
centers_5d <- matrix(rnorm(20), nrow = 4, ncol = 5)
r <- simulate_clusters(n = 400, centers = centers_5d, sds = 0.5, seed = 99)
```

---

### 5. simulate_prediction() --- Regression with Categorical Predictors

Generates regression data with continuous predictors, optional categorical predictors, and known ground-truth coefficients. Richer than `simulate_regression()` because it supports factor variables and computes population R-squared.

**Signature:**

```r
simulate_prediction(
  n,                          # Integer. Sample size. Min: 5.
  coefs,                      # Named numeric vector. May include "(Intercept)".
  cat_levels     = NULL,      # Named list of character vectors. Factor levels.
  cat_effects    = NULL,      # Named list of numeric vectors. Category effects.
                              #   NULL -> random effects from N(0, 4).
  error_sd       = 1,         # Numeric (> 0). Residual SD.
  predictor_means = NULL,     # Named numeric vector or NULL. Predictor means.
                              #   NULL -> all zeros.
  predictor_sds   = NULL,     # Named numeric vector or NULL. Predictor SDs.
                              #   NULL -> all ones.
  seed            = NULL       # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `y` | numeric | Response variable (always first column) |
| predictor names | numeric | Continuous predictors (names from `coefs`) |
| category names | factor | Categorical predictors (names from `cat_levels`) |

| `$params` fields | Type | Description |
|---|---|---|
| `coefs` | named numeric vector | Linear coefficients (including intercept) |
| `cat_effects` | list of numeric vectors | Per-level effects for each factor |
| `error_sd` | numeric | Residual standard deviation |
| `predictor_means` | named numeric vector | Continuous predictor means |
| `predictor_sds` | named numeric vector | Continuous predictor SDs |
| `r_squared` | numeric | Population R-squared = signal_var / (signal_var + error_sd^2) |

**Examples:**

```r
# Simple regression
r <- simulate_prediction(
  n = 200,
  coefs = c("(Intercept)" = 5, x1 = 2, x2 = -1),
  error_sd = 1, seed = 1
)
summary(lm(y ~ x1 + x2, data = r$data))
r$params$r_squared

# With categorical predictors
r <- simulate_prediction(
  n = 500,
  coefs = c("(Intercept)" = 10, age = 0.5, income = 0.001),
  cat_levels = list(
    treatment = c("control", "drug_a", "drug_b"),
    sex = c("M", "F")
  ),
  cat_effects = list(
    treatment = c(0, 3, 7),
    sex = c(0, -2)
  ),
  error_sd = 3,
  predictor_means = c(age = 45, income = 50000),
  predictor_sds = c(age = 12, income = 20000),
  seed = 42
)
summary(lm(y ~ ., data = r$data))

# High signal-to-noise (low error_sd -> high R^2)
r <- simulate_prediction(n = 100, coefs = c(x1 = 5), error_sd = 0.1, seed = 7)
r$params$r_squared  # very high
```

---

### 6. simulate_longitudinal() --- Panel / ESM / VAR(1) Data

Generates multilevel time-series data for testing longitudinal models: mlVAR, RI-CLPM, latent growth curves, dynamic structural equation models.

The data-generating process:
- **Temporal**: y(t) = mu_i + B * (y(t-1) - mu_i) + innovation
- **Contemporaneous**: innovations are correlated within time points
- **Between-person**: person-specific means mu_i ~ N(grand_means, between)

**Signature:**

```r
simulate_longitudinal(
  n               = 50L,          # Integer. Number of subjects. Min: 2.
  tp              = 50L,          # Integer. Time points per subject. Min: 3.
  vars            = 4L,           # Integer or character vector. Number/names of variables.
                                  #   Min: 2. Integer -> "V1", "V2", ...
  temporal        = NULL,         # p x p matrix or NULL. VAR(1) coefficient matrix B.
                                  #   B[i,j] = effect of var j at t-1 on var i at t.
                                  #   NULL -> auto-generated from ar_range, cross_range.
  contemporaneous = NULL,         # p x p matrix or NULL. Innovation correlation.
                                  #   NULL -> identity (uncorrelated).
  between         = NULL,         # p x p matrix or NULL. Person-mean covariance.
                                  #   NULL -> identity.
  grand_means     = NULL,         # Numeric vector[p] or NULL. Population means.
                                  #   NULL -> all zeros.
  innovation_sd   = 1,            # Numeric scalar or vector[p]. Innovation SDs.
  beeps_per_day   = NULL,         # Integer or NULL. ESM structure.
                                  #   tp must be divisible by beeps_per_day.
                                  #   NULL -> single "time" column.
  ar_range        = c(0.2, 0.5), # Numeric[2]. Autoregressive range (when temporal=NULL).
  cross_range     = c(0.05, 0.2),# Numeric[2]. Cross-lag range (when temporal=NULL).
  n_cross         = NULL,         # Integer or NULL. Number of cross-lags (when temporal=NULL).
                                  #   NULL -> min(3, p*(p-1)).
  complexity      = "clean",      # "clean", "auto", or character vector.
                                  #   Valid: "na", "outliers", "heavy_tailed",
                                  #   "heteroscedastic", "tiny_n".
  seed            = NULL           # Integer or NULL. Random seed.
                                  #   NULL -> auto-generated seed, stored in $seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `id` | integer | Subject identifier, 1 to n |
| `time` | integer | Time point, 1 to tp (when `beeps_per_day = NULL`) |
| `day` | integer | Day index (when `beeps_per_day` is set) |
| `beep` | integer | Beep within day (when `beeps_per_day` is set) |
| `V1`, ..., `Vp` | numeric | VAR(1) observed variables |

| `$params` fields | Type | Description |
|---|---|---|
| `temporal` | p x p matrix | VAR(1) coefficient matrix B (possibly rescaled) |
| `contemporaneous` | p x p matrix | Innovation correlation matrix |
| `between` | p x p matrix | Between-person covariance |
| `grand_means` | named numeric vector | Population means |
| `innovation_sd` | named numeric vector | Innovation standard deviations |
| `n` | integer | Number of subjects |
| `tp` | integer | Time points per subject |
| `var_names` | character vector | Variable names |

**Key behaviors:**
- Non-stationary B is auto-rescaled: `B <- B * (0.95 / max(Mod(eigen(B)$values)))`
- ESM mode: first beep of each day resets (no carry-over from previous day)
- `seed = NULL` auto-generates and stores a seed so results are always reproducible after the fact

**Examples:**

```r
# Basic: 50 subjects, 100 time points, 3 variables
r <- simulate_longitudinal(n = 50, tp = 100, vars = 3, seed = 42)
head(r$data)
r$params$temporal    # the true VAR(1) matrix

# Explicit temporal structure
B <- matrix(c(0.4,  0.0, 0.1,
              0.2,  0.3, 0.0,
              0.0, -0.1, 0.5), nrow = 3, byrow = TRUE)
r <- simulate_longitudinal(n = 80, tp = 100, vars = 3, temporal = B, seed = 1)

# Named variables
r <- simulate_longitudinal(
  n = 40, tp = 80,
  vars = c("mood", "energy", "stress"),
  grand_means = c(5, 6, 4),
  seed = 7
)

# ESM structure: 7 beeps/day, 10 days = 70 time points
r <- simulate_longitudinal(
  n = 50, tp = 70, vars = 4,
  beeps_per_day = 7, seed = 99
)
head(r$data)  # columns: id, day, beep, V1, V2, V3, V4

# Full explicit specification for parameter recovery
B <- matrix(c(0.3, 0.1, 0.1, 0.4), nrow = 2)
C <- matrix(c(1.0, 0.4, 0.4, 1.0), nrow = 2)
Bw <- matrix(c(2.0, 0.5, 0.5, 1.5), nrow = 2)
r <- simulate_longitudinal(
  n = 200, tp = 100, vars = 2,
  temporal = B,
  contemporaneous = C,
  between = Bw,
  grand_means = c(5, 3),
  innovation_sd = c(0.8, 1.2),
  seed = 42
)

# With edge-case injection
r <- simulate_longitudinal(
  n = 30, tp = 50, vars = 3,
  complexity = c("na", "outliers"),
  seed = 5
)
```

---

### 7. simulate_lpa() --- Latent Profile Analysis

Generates continuous-indicator data from a known latent profile (mixture of normals) structure.

**Signature:**

```r
simulate_lpa(
  means,          # Numeric matrix (n_vars x n_profiles). Column = profile means.
  sds,            # Scalar, vector, or matrix (n_vars x n_profiles). SDs. Recycled.
  props,          # Numeric vector[K]. Mixing proportions. Normalised internally.
  n,              # Integer. Total sample size.
  seed = NULL      # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `y1`, ..., `yp` | numeric | Observed indicator variables |
| `true_profile` | integer | True profile membership (1 to K) |

| `$params` fields | Type | Description |
|---|---|---|
| `means` | matrix (n_vars x n_profiles) | Input means |
| `sds` | matrix (n_vars x n_profiles) | Full SD matrix |
| `props` | numeric vector | Input proportions (un-normalised) |

**Examples:**

```r
# Two profiles, 4 indicators
means <- matrix(c(2, 2, 2, 2,     # profile 1: low on everything
                  8, 8, 8, 8),    # profile 2: high on everything
                nrow = 4, ncol = 2)
r <- simulate_lpa(means = means, sds = 1, props = c(0.6, 0.4), n = 500, seed = 1)

# Verify recovery: group means should match
tapply(r$data$y1, r$data$true_profile, mean)  # ~2 and ~8

# Three profiles with varying SDs
means3 <- matrix(c(0, 0,  5, 5,  10, 0), nrow = 2, ncol = 3)
sds3   <- matrix(c(0.5, 0.5,  1, 1,  2, 2), nrow = 2, ncol = 3)
r <- simulate_lpa(means = means3, sds = sds3, props = c(1, 2, 1), n = 1000, seed = 42)
```

---

### 8. simulate_lca() --- Latent Class Analysis

Generates binary-indicator data from a known latent class structure.

**Signature:**

```r
simulate_lca(
  item_probs,     # Numeric matrix (n_items x n_classes). P(item=1|class).
                  #   All values in [0, 1].
  class_probs,    # Numeric vector[K]. Class mixing proportions. Normalised internally.
  n,              # Integer. Total sample size.
  seed = NULL      # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `item1`, ..., `itemm` | integer (0/1) | Binary item responses |
| `true_class` | integer | True class membership (1 to K) |

| `$params` fields | Type | Description |
|---|---|---|
| `item_probs` | matrix (n_items x n_classes) | Input response probabilities |
| `class_probs` | numeric vector | Input proportions (un-normalised) |

**Examples:**

```r
# Two classes, 6 items: class 1 endorses items 1-3, class 2 endorses items 4-6
item_probs <- matrix(c(
  0.9, 0.1,
  0.8, 0.2,
  0.9, 0.1,
  0.1, 0.9,
  0.2, 0.8,
  0.1, 0.9
), nrow = 6, ncol = 2, byrow = TRUE)
r <- simulate_lca(item_probs = item_probs, class_probs = c(0.5, 0.5),
                   n = 500, seed = 1)
colMeans(r$data[r$data$true_class == 1, 1:6])  # high on items 1-3
colMeans(r$data[r$data$true_class == 2, 1:6])  # high on items 4-6
```

---

### 9. simulate_regression() --- Simple Linear Regression

Generates regression data with known coefficients and independent normal predictors. For regression with categorical predictors, use `simulate_prediction()` instead.

**Signature:**

```r
simulate_regression(
  coefs,            # Named numeric vector. Must include names. May include "(Intercept)".
  predictor_sds,    # Named numeric vector. Names must match non-intercept coefs. All > 0.
  error_sd,         # Numeric (> 0). Residual SD.
  n,                # Integer. Sample size.
  seed = NULL        # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `y` | numeric | Response (first column) |
| predictor names | numeric | Independent predictors |

| `$params` fields | Type | Description |
|---|---|---|
| `coefs` | named numeric vector | Input coefficients |
| `predictor_sds` | named numeric vector | Input predictor SDs |
| `error_sd` | numeric | Input residual SD |

**Examples:**

```r
coefs <- c("(Intercept)" = 2, x1 = 3, x2 = -1)
r <- simulate_regression(
  coefs = coefs,
  predictor_sds = c(x1 = 1, x2 = 1),
  error_sd = 0.5, n = 500, seed = 42
)
coef(lm(y ~ ., data = r$data))  # should be close to c(2, 3, -1)
```

---

### 10. simulate_fa() --- Factor Analysis

Generates multivariate normal data from an explicit factor model: Sigma = Lambda * Phi * Lambda' + Psi.

**Signature:**

```r
simulate_fa(
  loadings,        # Numeric matrix (p x m). Factor loadings. Must be matrix.
  phi = NULL,      # Numeric matrix (m x m) or NULL. Factor correlations.
                   #   NULL -> diag(m) (orthogonal factors).
  psi = NULL,      # Numeric vector[p] or NULL. Unique variances.
                   #   NULL -> auto-computed as 1 - communalities.
  n,               # Integer. Sample size.
  seed = NULL       # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `y1`, ..., `yp` | numeric | Observed variables |

| `$params` fields | Type | Description |
|---|---|---|
| `loadings` | matrix (p x m) | Factor loadings |
| `phi` | matrix (m x m) | Factor correlation matrix |
| `psi` | numeric vector | Unique variances |
| `sigma_implied` | matrix (p x p) | Model-implied covariance: Lambda * Phi * Lambda' + Psi |

**Examples:**

```r
# Two factors, 6 items (3 per factor), orthogonal
loadings <- matrix(c(
  0.8, 0.0,
  0.7, 0.0,
  0.6, 0.0,
  0.0, 0.8,
  0.0, 0.7,
  0.0, 0.6
), nrow = 6, ncol = 2, byrow = TRUE)
r <- simulate_fa(loadings = loadings, n = 500, seed = 1)
factanal(r$data, factors = 2)

# Oblique (correlated factors)
phi <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
r <- simulate_fa(loadings = loadings, phi = phi, n = 500, seed = 42)

# Compare implied vs observed covariance
r$params$sigma_implied
cov(r$data)
```

---

### 11. simulate_seq_clusters() --- Markov Chain Sequence Clusters

Generates wide-format categorical sequence data where each row is drawn from one of K Markov chains with distinct transition matrices.

**Signature:**

```r
simulate_seq_clusters(
  trans_list  = NULL,      # List of square row-stochastic matrices, or NULL.
  props       = NULL,      # Numeric vector[K] or NULL. Mixing proportions.
  n           = 300L,      # Integer. Number of sequences.
  seq_length  = 20L,       # Integer. Time points per sequence.
  init_probs  = NULL,      # Vector, list of vectors, or NULL. Initial state probs.
  n_clusters  = 3L,        # Integer. Clusters when trans_list=NULL.
  n_states    = 10L,       # Integer. States when trans_list=NULL.
  states      = NULL,      # Character vector or NULL. State names.
  seed        = NULL        # Integer or NULL. Random seed.
)
```

**Returns `saqr_sim` with:**

| `$data` columns | Type | Description |
|---|---|---|
| `T1`, ..., `T{seq_length}` | character | State labels at each time point |
| `true_cluster` | integer | True cluster membership (1 to K) |

| `$params` fields | Type | Description |
|---|---|---|
| `trans_list` | list of matrices | Per-cluster transition matrices |
| `props` | numeric vector | Normalised mixing proportions |
| `init_probs` | list of numeric vectors | Per-cluster initial probabilities |

**Examples:**

```r
# Auto-generated: 3 clusters, 10 states
r <- simulate_seq_clusters(n = 500, n_clusters = 3, n_states = 5, seed = 1)
table(r$data$true_cluster)

# Explicit transition matrices
m1 <- matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2,
             dimnames = list(c("A","B"), c("A","B")))
m2 <- matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2,
             dimnames = list(c("A","B"), c("A","B")))
r <- simulate_seq_clusters(trans_list = list(m1, m2),
                            props = c(0.6, 0.4), n = 200, seed = 42)
```

---

### 12. simulate_data() --- Random-Parameter Quick Generator

Generates datasets with **randomised** structural parameters (sample sizes, effect sizes, group counts). The seed determines both the random data and the structural parameters. Use this for stress testing and robustness checks, not parameter recovery.

**Signature:**

```r
simulate_data(
  type,                       # Character. One of: "ttest", "anova", "correlation",
                              #   "clusters", "factor_analysis", "prediction",
                              #   "mlvar", "batch".
  seed       = NULL,          # Integer or NULL. Determines structure + data.
  complexity = "clean",       # "clean", "auto", or character vector of edge cases.
  ...,                        # Named overrides (n, n_groups, effect_size, etc.).
  n_batch    = NULL            # Integer or NULL. Batch mode.
)
```

**Returns:**

- **Single dataset**: `data.frame` with attributes `type`, `info` (suggested analysis code), `complexity`.
- **Batch mode**: list of data.frames (class `sim_batch_type` or `sim_batch_full`).

**Note:** This function does NOT return `saqr_sim`. It returns bare data.frames because the parameters are random (not user-specified). For explicit-parameter simulation, use the dedicated functions above.

**Complexity options:** `"na"`, `"outliers"`, `"ties"`, `"duplicates"`, `"constant_col"`, `"all_na_col"`, `"tiny_n"`, `"heavy_tailed"`, `"heteroscedastic"`, `"extreme_imbalance"`, `"multicollinear"`.

**Type output structure:**

| Type | Columns | Attributes |
|---|---|---|
| `"ttest"` | `group` (factor), `score` | --- |
| `"anova"` | `group` (factor), `score` | --- |
| `"correlation"` | `x1`--`xp` (4--7 numeric) | --- |
| `"clusters"` | `x1`--`xd`, `true_cluster` | --- |
| `"factor_analysis"` | `x1`--`xp` (6--16 numeric) | `n_factors`, `loadings` |
| `"prediction"` | `y`, `x1`--`x4`, `cat1`, `cat2` | --- |
| `"mlvar"` | `id`, `day`, `beep`, `V1`--`Vd` | `true_temporal`, `true_contemporaneous`, `vars` |

**Examples:**

```r
# Quick dataset for any analysis
d <- simulate_data("ttest", seed = 42)
t.test(score ~ group, data = d)

# With edge-case injection
d <- simulate_data("correlation", seed = 1, complexity = c("na", "outliers"))

# Batch: 100 datasets for simulation studies
batch <- simulate_data("ttest", seed = 1, n_batch = 100)
p_values <- vapply(batch, function(d) {
  t.test(score ~ group, data = d)$p.value
}, numeric(1))

# All 7 types at once
all_data <- simulate_data("batch", seed = 1, n_batch = 50)
```

---

## Parameter Recovery Cookbook

### t-test: Recover means and Cohen's d

```r
r <- simulate_ttest(n_a = 500, n_b = 500, mean_a = 50, mean_b = 55,
                     sd_a = 10, sd_b = 10, seed = 42)
fit <- t.test(score ~ group, data = r$data)

# Compare:
true_diff   <- r$params$mean_b - r$params$mean_a   # 5
obs_diff    <- diff(fit$estimate)                   # should be ~5
true_d      <- r$params$cohens_d                    # 0.5
```

### ANOVA: Recover group means and eta-squared

```r
r <- simulate_anova(n = 200, means = c(10, 15, 20), sds = 3, seed = 1)
fit <- aov(score ~ group, data = r$data)

# Compare:
true_means  <- r$params$means
obs_means   <- tapply(r$data$score, r$data$group, mean)
true_eta    <- r$params$eta_squared
obs_eta     <- summary(fit)[[1]]["group", "Sum Sq"] /
               sum(summary(fit)[[1]][, "Sum Sq"])
```

### Correlation: Recover correlation matrix

```r
R <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.4, 0.2, 0.4, 1), nrow = 3)
r <- simulate_correlation(n = 1000, sigma = R, seed = 42)

# Compare:
true_cor <- r$params$sigma
obs_cor  <- cor(r$data)
max(abs(true_cor - obs_cor))  # should be small
```

### Regression: Recover coefficients and R-squared

```r
r <- simulate_prediction(
  n = 2000,
  coefs = c("(Intercept)" = 10, x1 = 3, x2 = -2),
  error_sd = 1, seed = 1
)
fit <- lm(y ~ x1 + x2, data = r$data)

# Compare:
true_coefs <- r$params$coefs
obs_coefs  <- coef(fit)
true_r2    <- r$params$r_squared
obs_r2     <- summary(fit)$r.squared
```

### Factor Analysis: Recover loadings

```r
L <- matrix(c(0.8, 0.7, 0.6, 0, 0, 0,
              0, 0, 0, 0.8, 0.7, 0.6), nrow = 6, ncol = 2, byrow = FALSE)
r <- simulate_fa(loadings = L, n = 2000, seed = 42)
fit <- factanal(r$data, factors = 2, rotation = "varimax")

# Compare:
true_loadings <- r$params$loadings
obs_loadings  <- fit$loadings[]
# Note: rotation may permute/reflect factors
```

### Clustering: Recover centers

```r
true_centers <- matrix(c(0, 0, 5, 5, 10, 0), nrow = 3, byrow = TRUE)
r <- simulate_clusters(n = 3000, centers = true_centers, sds = 0.5, seed = 1)
fit <- kmeans(r$data[, c("x1", "x2")], centers = 3, nstart = 25)

# Compare (after matching clusters):
obs_centers <- fit$centers[order(fit$centers[, 1]), ]
true_sorted <- true_centers[order(true_centers[, 1]), ]
max(abs(obs_centers - true_sorted))  # should be small
```

### mlVAR: Recover temporal matrix

```r
B <- matrix(c(0.4, 0.0,
              0.2, 0.3), nrow = 2, byrow = TRUE)
r <- simulate_longitudinal(
  n = 200, tp = 100, vars = 2,
  temporal = B,
  between = diag(2) * 0.5,
  innovation_sd = 0.5, seed = 42
)

# Person-mean center and compute lagged variables
d <- r$data
d$V1c <- d$V1 - ave(d$V1, d$id, FUN = mean)
d$V2c <- d$V2 - ave(d$V2, d$id, FUN = mean)
d$V1_lag <- ave(d$V1c, d$id, FUN = function(x) c(NA, x[-length(x)]))
d$V2_lag <- ave(d$V2c, d$id, FUN = function(x) c(NA, x[-length(x)]))
d <- d[complete.cases(d), ]

# Row 1 of B: V1 ~ V1_lag + V2_lag
fit1 <- lm(V1c ~ V1_lag + V2_lag - 1, data = d)
coef(fit1)  # should approximate c(0.4, 0.0)

# Row 2 of B: V2 ~ V1_lag + V2_lag
fit2 <- lm(V2c ~ V1_lag + V2_lag - 1, data = d)
coef(fit2)  # should approximate c(0.2, 0.3)
```

### LPA: Recover profile means

```r
means <- matrix(c(2, 2, 8, 8), nrow = 2, ncol = 2)
r <- simulate_lpa(means = means, sds = 0.5, props = c(0.5, 0.5), n = 2000, seed = 1)

# Compare per-profile means:
true_means <- r$params$means
obs_means  <- vapply(1:2, function(k) {
  colMeans(r$data[r$data$true_profile == k, c("y1", "y2")])
}, numeric(2))
max(abs(true_means - obs_means))  # should be small
```

---

## Summary Table

| Function | Use case | Key params | Computed ground truth |
|---|---|---|---|
| `simulate_ttest` | Two-group comparison | `mean_a`, `mean_b`, `sd_a`, `sd_b` | `cohens_d` |
| `simulate_anova` | Multi-group comparison | `means`, `sds`, `n` | `eta_squared` |
| `simulate_correlation` | Multivariate association | `sigma` | `is_correlation` |
| `simulate_clusters` | Cluster analysis | `centers`, `sds`, `props` | `true_cluster` in data |
| `simulate_prediction` | Regression + categories | `coefs`, `cat_levels`, `cat_effects` | `r_squared` |
| `simulate_longitudinal` | mlVAR / ESM / panel | `temporal`, `contemporaneous`, `between` | all three matrices |
| `simulate_lpa` | Latent Profile Analysis | `means`, `sds`, `props` | `true_profile` in data |
| `simulate_lca` | Latent Class Analysis | `item_probs`, `class_probs` | `true_class` in data |
| `simulate_regression` | Simple linear regression | `coefs`, `predictor_sds`, `error_sd` | exact inputs |
| `simulate_fa` | Factor analysis | `loadings`, `phi`, `psi` | `sigma_implied` |
| `simulate_seq_clusters` | Sequence clustering | `trans_list`, `props` | `true_cluster` in data |
| `simulate_data` | Quick/stress testing | `type`, `complexity` | randomised (in attrs) |

All functions except `simulate_data()` return `saqr_sim` objects. Access `$data` for the data.frame, `$params` for ground-truth parameters, `$seed` for reproducibility.
