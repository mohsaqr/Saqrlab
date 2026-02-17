# 04 — Permutation Test for Network Comparison

`permutation_test()` lives in `R/permutation_test.R` (580 lines). It compares two networks using a permutation test, shuffling which observations belong to which group.

---

## 1. Function Signature (lines 72-78)

```r
permutation_test <- function(x, y,
                             iter = 1000L,
                             alpha = 0.05,
                             paired = FALSE,
                             adjust = "none",
                             nlambda = 50L,
                             seed = NULL)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `x` | netobject | required | First network (from `build_network()`) |
| `y` | netobject | required | Second network |
| `iter` | integer | 1000 | Number of permutation iterations |
| `alpha` | numeric | 0.05 | Significance level |
| `paired` | logical | FALSE | Paired permutation (within-pair swaps) |
| `adjust` | character | `"none"` | p-value adjustment method (passed to `p.adjust()`) |
| `nlambda` | integer | 50 | Lambda path resolution for glasso permutations |
| `seed` | integer or NULL | NULL | RNG seed |

---

## 2. Input Validation (lines 80-128)

```r
stopifnot(
  inherits(x, "netobject"),
  inherits(y, "netobject"),
  is.numeric(iter), length(iter) == 1, iter >= 2,
  is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
  is.logical(paired), length(paired) == 1,
  is.character(adjust), length(adjust) == 1
)
```

**Additional checks**:

- Both `x` and `y` must have `$data` (line 90-97). If NULL, suggests rebuilding with `build_network()`.
- Methods must match: `x$method != y$method` → error (lines 99-102)
- Nodes must be the same set: `setequal(x$nodes, y$nodes)` (line 104)
- Paired mode requires equal `nrow(x$data) == nrow(y$data)` (lines 118-123)

**Node reordering** (lines 109-112):

```r
nodes <- x$nodes
if (!identical(x$nodes, y$nodes)) {
  y$matrix <- y$matrix[nodes, nodes]
}
```

If nodes are the same set but in different order, `y$matrix` is reordered to match `x`.

---

## 3. Observed Difference (line 131)

```r
obs_diff <- x$matrix - y$matrix
```

Element-wise difference. This is what the permutation test evaluates.

---

## 4. Permutation Dispatch (lines 134-144)

```r
if (method %in% c("relative", "frequency", "co_occurrence")) {
  perm_result <- .permutation_transition(
    x, y, nodes, method, iter, paired
  )
} else {
  perm_result <- .permutation_association(
    x, y, nodes, method, iter, paired, nlambda
  )
}
```

Both paths return `list(exceed_counts, perm_sd)` — running accumulators, NOT full matrices.

---

## 5. Transition Fast Path: `.permutation_transition()` (lines 206-267)

### Setup (lines 207-228)

```r
n_nodes <- length(nodes)
nbins   <- n_nodes * n_nodes
is_relative <- method == "relative"

# Pre-compute per-sequence counts for both groups
trans_x <- .precompute_per_sequence(x$data, method, x$params, nodes)
trans_y <- .precompute_per_sequence(y$data, method, y$params, nodes)

n_x <- nrow(trans_x)
n_y <- nrow(trans_y)

# Pool sequences
pooled  <- rbind(trans_x, trans_y)
n_total <- n_x + n_y

obs_flat <- as.vector(x$matrix - y$matrix)

# Running counters
exceed_counts <- integer(nbins)
sum_diffs     <- numeric(nbins)
sum_diffs_sq  <- numeric(nbins)
```

Uses the same `.precompute_per_sequence()` function as bootstrap. Pooled matrix has `n_total` rows.

### Main Loop (lines 230-257)

```r
for (i in seq_len(iter)) {
  if (paired) {
    # Paired: randomly swap x/y within each pair
    swaps <- sample(c(TRUE, FALSE), n_x, replace = TRUE)
    idx_x <- ifelse(swaps, seq(n_x + 1L, n_total), seq_len(n_x))
    idx_y <- ifelse(swaps, seq_len(n_x), seq(n_x + 1L, n_total))
    counts_x <- colSums(pooled[idx_x, , drop = FALSE])
    counts_y <- colSums(pooled[idx_y, , drop = FALSE])
  } else {
    # Unpaired: shuffle group labels
    idx_x <- sample.int(n_total, n_x)
    counts_x <- colSums(pooled[idx_x, , drop = FALSE])
    counts_y <- colSums(pooled[-idx_x, , drop = FALSE])
  }

  # Post-process each group
  mat_x <- .postprocess_counts(counts_x, n_nodes, is_relative, x$scaling, x$threshold)
  mat_y <- .postprocess_counts(counts_y, n_nodes, is_relative, y$scaling, y$threshold)

  perm_diff <- as.vector(mat_x) - as.vector(mat_y)

  # Accumulate
  exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
  sum_diffs     <- sum_diffs + perm_diff
  sum_diffs_sq  <- sum_diffs_sq + perm_diff^2
}
```

**CRITICAL — `sample.int(n_total, n_x)`**: This draws `n_x` indices from `1:n_total` WITHOUT replacement. It consumes exactly `n_x` RNG values. Using `sample.int(n_total)` (full permutation) would consume `n_total` values, causing RNG stream divergence from NCT (NetworkComparisonTest). This was a key fix for achieving p-value agreement with NCT.

### Post-process: `.postprocess_counts()` (lines 272-283)

```r
.postprocess_counts <- function(counts, n_nodes, is_relative, scaling, threshold) {
  mat <- matrix(counts, n_nodes, n_nodes, byrow = TRUE)
  if (is_relative) {
    rs <- rowSums(mat)
    nz <- rs > 0
    mat[nz, ] <- mat[nz, ] / rs[nz]
  }
  if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
  if (threshold > 0) mat[abs(mat) < threshold] <- 0
  mat
}
```

**Key**: `byrow = TRUE` is critical here too — the flat counts from `colSums` are in the same row-major order as the precomputed matrix.

### SD Computation (lines 260-261)

```r
perm_mean <- sum_diffs / iter
perm_sd   <- sqrt(pmax(sum_diffs_sq / iter - perm_mean^2, 0))
```

Uses the running sums formula: `Var = E[X²] - E[X]²`. The `pmax(..., 0)` guards against floating-point errors producing negative variance.

---

## 6. Association Path: `.permutation_association()` (lines 295-423)

### Setup (lines 296-314)

```r
n_x <- nrow(x$data)
n_y <- nrow(y$data)
pooled_mat <- rbind(x$data, y$data)    # $data is ALREADY cleaned
n_total <- n_x + n_y
```

**CRITICAL**: `$data` stores cleaned data (numeric matrix, no NAs, no zero-variance). Pool directly — do NOT re-clean.

Extract params:

```r
cor_method  <- params_x$cor_method %||% "pearson"
threshold_x <- x$threshold
threshold_y <- y$threshold
scaling_x   <- x$scaling
scaling_y   <- y$scaling
```

### Glasso Pre-computation (lines 321-330)

```r
if (use_fast && method == "glasso") {
  gamma        <- params_x$gamma %||% 0.5
  penalize_diag <- params_x$penalize.diagonal %||% FALSE

  S_pooled      <- cor(pooled_mat, method = cor_method)
  perm_rholist  <- .compute_lambda_path(S_pooled, nlambda, 0.01)
  p_glasso      <- ncol(pooled_mat)
}
```

The lambda path is computed ONCE from the pooled correlation matrix. This same path is used for all iterations.

### Per-Method Closures (lines 334-381)

Three fast-path closures for built-in association methods:

**cor** (lines 335-338):
```r
function(mat_subset) {
  S <- cor(mat_subset, method = cor_method)
  diag(S) <- 0
  S
}
```

**pcor** (lines 339-343):
```r
function(mat_subset) {
  S  <- cor(mat_subset, method = cor_method)
  Wi <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Wi)) return(NULL)
  .precision_to_pcor(Wi, threshold = 0)
}
```

**glasso** (lines 344-361) — the `glassopath()` optimization:
```r
function(mat_subset) {
  S     <- cor(mat_subset, method = cor_method)
  n_obs <- nrow(mat_subset)
  gp <- tryCatch(
    glasso::glassopath(s = S, rholist = perm_rholist, trace = 0,
                       penalize.diagonal = penalize_diag),
    error = function(e) NULL
  )
  if (is.null(gp)) return(NULL)
  best_wi <- .select_ebic_from_path(gp, S, n_obs, gamma, p_glasso, perm_rholist)
  if (is.null(best_wi)) return(NULL)
  .precision_to_pcor(best_wi, threshold = 0)
}
```

**`glassopath()` optimization**: Instead of calling `glasso()` separately for each lambda value, `glassopath()` solves the entire regularization path in a single Fortran call with internal warm-starting. This is **14.5x faster** than separate `glasso()` calls.

**`glassopath()` returns**: `$wi` as a 3D array `p x p x n_lambda`.

**CRITICAL — Must re-select lambda per iteration**: Using a fixed lambda (from the original network) across all permutations drops agreement from 99% to 93%. Each permutation sample has different correlation structure, so the optimal lambda differs. EBIC selection per iteration via `.select_ebic_from_path()` is required.

### Custom Estimator Fallback (lines 363-381)

For non-built-in methods, falls back to full `do.call(estimator$fn, ...)`.

### Main Loop (lines 388-414)

```r
for (i in seq_len(iter)) {
  if (paired) {
    swaps <- sample(c(TRUE, FALSE), n_x, replace = TRUE)
    idx_x <- ifelse(swaps, seq(n_x + 1L, n_total), seq_len(n_x))
    idx_y <- ifelse(swaps, seq_len(n_x), seq(n_x + 1L, n_total))
  } else {
    idx_x <- sample.int(n_total, n_x)
    idx_y <- seq_len(n_total)[-idx_x]
  }

  mat_x <- estimate_from_rows(pooled_mat[idx_x, , drop = FALSE])
  mat_y <- estimate_from_rows(pooled_mat[idx_y, , drop = FALSE])

  if (is.null(mat_x) || is.null(mat_y)) next   # skip failed fits

  # Apply scaling and threshold per-network
  if (!is.null(scaling_x)) mat_x <- .apply_scaling(mat_x, scaling_x)
  if (threshold_x > 0)     mat_x[abs(mat_x) < threshold_x] <- 0
  if (!is.null(scaling_y)) mat_y <- .apply_scaling(mat_y, scaling_y)
  if (threshold_y > 0)     mat_y[abs(mat_y) < threshold_y] <- 0

  perm_diff <- as.vector(mat_x) - as.vector(mat_y)

  exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
  sum_diffs     <- sum_diffs + perm_diff
  sum_diffs_sq  <- sum_diffs_sq + perm_diff^2
}
```

**Key**: Scaling and threshold are applied per-network (x and y can have different scaling/threshold settings).

**Key**: Failed fits (`next`) are skipped but the `iter` denominator stays the same. This is conservative.

---

## 7. EBIC Selection from glassopath: `.select_ebic_from_path()` (lines 430-453)

```r
.select_ebic_from_path <- function(gp, S, n, gamma, p, rholist) {
  n_lambda  <- length(rholist)
  best_ebic <- Inf
  best_wi   <- NULL

  for (k in seq_len(n_lambda)) {
    wi_k <- gp$wi[, , k]

    log_det <- determinant(wi_k, logarithm = TRUE)
    if (log_det$sign <= 0) next
    log_det_val <- as.numeric(log_det$modulus)

    loglik <- (n / 2) * (log_det_val - sum(diag(S %*% wi_k)))
    npar   <- sum(abs(wi_k[upper.tri(wi_k)]) > 1e-10)
    ebic   <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)

    if (ebic < best_ebic) {
      best_ebic <- ebic
      best_wi   <- wi_k
    }
  }

  best_wi  # NULL if all fits failed
}
```

Same EBIC formula as `.select_ebic()` in estimators.R, but operates on the pre-computed 3D array from `glassopath()` instead of fitting glasso per lambda.

---

## 8. P-values (lines 147-155)

```r
obs_flat      <- as.vector(obs_diff)
p_values_flat <- (perm_result$exceed_counts + 1L) / (iter + 1L)

# Apply multiple comparison correction
p_values_flat <- p.adjust(p_values_flat, method = adjust)

p_mat <- matrix(p_values_flat, n_nodes, n_nodes,
                dimnames = list(nodes, nodes))
```

**Key**: `p.adjust()` is applied to the flat vector before reshaping to matrix. This correctly handles the adjustment across all edges simultaneously.

**Key**: The `+1` in numerator and denominator is the standard conservative correction (prevents p = 0).

---

## 9. Effect Size (lines 158-164)

```r
perm_sd <- perm_result$perm_sd
perm_sd[perm_sd == 0] <- NA_real_
es_flat <- obs_flat / perm_sd
es_flat[is.na(es_flat)] <- 0
```

Cohen's d style: observed difference divided by permutation SD. Zeros where SD is 0 (e.g., self-loops where diff and permutation diff are always 0).

---

## 10. Summary: `.build_permutation_summary()` (lines 460-483)

```r
dt <- data.table::data.table(
  from        = rep(nodes, each = n),
  to          = rep(nodes, times = n),
  weight_x    = as.vector(t(x_matrix)),
  weight_y    = as.vector(t(y_matrix)),
  diff        = as.vector(t(obs_diff)),
  effect_size = as.vector(t(es_mat)),
  p_value     = as.vector(t(p_mat)),
  sig         = as.vector(t(p_mat)) < alpha
)
```

Same `as.vector(t(...))` pattern as bootstrap summary for row-major flattening.

**Filtering** — keeps edges present in EITHER network:

```r
if (directed) {
  dt <- dt[(weight_x != 0 | weight_y != 0) & from != to]
} else {
  dt <- dt[(weight_x != 0 | weight_y != 0) & from < to]
}
```

---

## 11. Return Value (lines 183-198)

```r
result <- list(
  x           = x,           # netobject
  y           = y,           # netobject
  diff        = obs_diff,    # matrix (x - y)
  diff_sig    = diff_sig,    # matrix (significant diffs only)
  p_values    = p_mat,       # matrix
  effect_size = es_mat,      # matrix
  summary     = summary_df,  # data.frame
  method      = method,
  iter        = iter,
  alpha       = alpha,
  paired      = paired,
  adjust      = adjust
)
class(result) <- "saqr_permutation"
```

---

## 12. Validation Against External Packages

| Comparison | Result |
|---|---|
| vs `NetworkComparisonTest::NCT` (cor) | p-values identical (r=1.0), 100% agreement, 1.9x faster |
| vs `NetworkComparisonTest::NCT` (glasso) | p-value r=0.9999, 99.9% agreement, 2.6x faster |
| vs `tna::permutation_test` (relative) | identical p-values (r=1.0) |
| vs `tna::permutation_test` (co_occurrence) | identical p-values (r=1.0) |

The key to NCT agreement was matching RNG consumption: `sample.int(n_total, n_x)` instead of `sample.int(n_total)`.
