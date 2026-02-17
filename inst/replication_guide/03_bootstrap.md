# 03 — Bootstrap Network Estimation

`bootstrap_network()` lives in `R/bootstrap_network.R` (599 lines). It provides non-parametric bootstrap for any network estimated by `build_network()`.

---

## 1. Function Signature (lines 84-96)

```r
bootstrap_network <- function(data,
                              method = "relative",
                              params = list(),
                              scaling = NULL,
                              threshold = 0,
                              level = NULL,
                              id_col = NULL,
                              iter = 1000L,
                              ci_level = 0.05,
                              inference = "stability",
                              consistency_range = c(0.75, 1.25),
                              edge_threshold = NULL,
                              seed = NULL)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | data.frame or matrix | required | Same as `build_network()` |
| `method` | character | `"relative"` | Estimator name/alias |
| `params` | list | `list()` | Method-specific params |
| `scaling`, `threshold`, `level`, `id_col` | various | — | Passed to `build_network()` |
| `iter` | integer | 1000 | Number of bootstrap iterations |
| `ci_level` | numeric | 0.05 | Significance level for CIs and p-values |
| `inference` | character | `"stability"` | `"stability"` or `"threshold"` |
| `consistency_range` | numeric(2) | `c(0.75, 1.25)` | Multiplicative bounds for stability inference |
| `edge_threshold` | numeric or NULL | NULL | Fixed threshold for threshold inference (auto-set if NULL) |
| `seed` | integer or NULL | NULL | RNG seed for reproducibility |

---

## 2. Input Validation (lines 97-113)

```r
stopifnot(
  is.character(method), length(method) == 1,
  is.list(params),
  is.numeric(threshold), length(threshold) == 1, threshold >= 0,
  is.numeric(iter), length(iter) == 1, iter >= 2,
  is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
  is.character(inference), length(inference) == 1,
  is.numeric(consistency_range), length(consistency_range) == 2
)
iter <- as.integer(iter)
inference <- match.arg(inference, c("stability", "threshold"))
```

**Key**: `iter >= 2` (not >= 1). Seed is set before any estimation work:

```r
if (!is.null(seed)) {
  stopifnot(is.numeric(seed), length(seed) == 1)
  set.seed(seed)
}
```

---

## 3. Compute Original Network (lines 116-128)

```r
method <- .resolve_method_alias(method)
estimator <- get_estimator(method)
directed <- estimator$directed

original <- build_network(
  data = data, method = method, params = params,
  scaling = scaling, threshold = threshold,
  level = level, id_col = id_col
)
```

The original network is computed first via `build_network()`. This ensures the original and bootstrap use identical estimation logic.

---

## 4. Bootstrap Dispatch (lines 131-142)

```r
if (method %in% c("relative", "frequency", "co_occurrence")) {
  boot_matrices <- .bootstrap_transition(
    data, method, params, states, scaling, threshold, iter
  )
} else {
  boot_matrices <- .bootstrap_association(
    data, estimator, params, states, scaling, threshold, iter,
    level, id_col
  )
}
```

Both paths return an `iter x nbins` matrix where `nbins = n_states^2`. Each row is a flattened network matrix from one bootstrap iteration.

---

## 5. Transition Fast Path (lines 227-357)

### 5.1 Top-level: `.bootstrap_transition()` (lines 227-254)

```r
.bootstrap_transition <- function(data, method, params, states,
                                  scaling, threshold, iter) {
  n_states <- length(states)
  nbins <- n_states * n_states
  is_relative <- method == "relative"

  # Pre-compute per-sequence count matrix ONCE
  trans_2d <- .precompute_per_sequence(data, method, params, states)
  n_seq <- nrow(trans_2d)

  # vapply: each iteration resamples sequences + sums + post-processes
  boot_flat <- vapply(seq_len(iter), function(i) {
    idx <- sample.int(n_seq, n_seq, replace = TRUE)
    boot_counts <- colSums(trans_2d[idx, , drop = FALSE])
    # Inline post-processing
    mat <- matrix(boot_counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    as.vector(mat)
  }, numeric(nbins))

  t(boot_flat)  # transpose: vapply returns nbins x iter → t() gives iter x nbins
}
```

**Key insight**: The entire bootstrap is done without calling `build_network()` or any estimator function per iteration. Per-sequence counts are precomputed once, then each iteration just resamples row indices and sums.

**Key**: `vapply()` returns an `nbins x iter` matrix (each call returns a vector of length `nbins`). `t()` transposes to `iter x nbins`.

**Key**: Post-processing is inlined (no dimnames on the matrix since it's flattened immediately). This avoids overhead of matrix naming per iteration.

### 5.2 Precomputation: `.precompute_per_sequence()` (lines 266-286)

Determines format (auto → wide/long). **Errors if long format** — bootstrap fast path requires wide data:

```r
if (format == "long") {
  stop("Bootstrap fast path requires wide-format data. ",
       "Convert with convert_sequence_format() first.", call. = FALSE)
}
```

Delegates to `.precompute_per_sequence_wide()`.

### 5.3 Per-Sequence Counts: `.precompute_per_sequence_wide()` (lines 294-357)

Returns an `n_seq x nbins` matrix where each row holds the flattened counts for one sequence.

**For transition methods (relative/frequency)** — lines 305-323:

```r
# Integer-encode the entire matrix once
int_mat <- matrix(match(mat, states), nrow = nr, ncol = nc)

# Consecutive column pairs
from_mat <- int_mat[, -nc, drop = FALSE]
to_mat   <- int_mat[, -1L, drop = FALSE]

row_ids  <- rep(seq_len(nr), times = nc - 1L)
from_vec <- as.vector(from_mat)
to_vec   <- as.vector(to_mat)

valid <- !is.na(from_vec) & !is.na(to_vec)
row_ids  <- row_ids[valid]
from_vec <- from_vec[valid]
to_vec   <- to_vec[valid]

pair_idx     <- (from_vec - 1L) * n_states + to_vec
combined_idx <- (row_ids - 1L) * nbins + pair_idx

counts <- tabulate(combined_idx, nbins = nr * nbins)
matrix(as.numeric(counts), nrow = nr, ncol = nbins, byrow = TRUE)
```

**Algorithm explained**: The combined flat index `(row_id - 1) * nbins + pair_idx` maps each transition pair to a unique position in an `nr x nbins` flat array. A single `tabulate()` call counts all pairs across all sequences simultaneously. Reshaping with `byrow = TRUE` gives per-sequence counts.

**For co-occurrence method** — lines 326-356:

Nested loop over column pairs `(i, j)` where `i < j`, accumulating bidirectional counts per sequence. Uses the same combined index pattern but with forward + reverse directions:

```r
combined <- c(
  (rows_valid - 1L) * nbins + idx_fwd,
  (rows_valid - 1L) * nbins + idx_rev
)
```

Self-pair halving applied at the end: `result[, diag_indices] <- result[, diag_indices] / 2`.

---

## 6. Association Path (lines 364-403)

### `.bootstrap_association()`

Full estimator call per iteration. No precomputation shortcut.

```r
boot_flat <- vapply(seq_len(iter), function(i) {
  boot_data <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]

  # For multilevel: apply decomposition to bootstrap sample
  if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
    boot_data <- tryCatch(
      .decompose_multilevel(boot_data, id_col = id_col, level = level),
      error = function(e) NULL
    )
    if (is.null(boot_data)) return(rep(NA_real_, nbins))
  }

  est <- tryCatch(
    do.call(estimator$fn, c(list(data = boot_data), params)),
    error = function(e) NULL
  )
  if (is.null(est)) return(rep(NA_real_, nbins))

  mat <- est$matrix
  # Align to expected states order
  if (!identical(rownames(mat), states)) {
    common <- intersect(states, rownames(mat))
    if (length(common) < n_states) return(rep(NA_real_, nbins))
    mat <- mat[states, states]
  }

  if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
  if (threshold > 0) mat[abs(mat) < threshold] <- 0
  as.vector(mat)
}, numeric(nbins))

t(boot_flat)
```

**Key**: Failed estimator fits return `NA_real_` vectors. These are handled in the statistics computation (filtered by `valid_rows`).

**Key**: State alignment — if the bootstrap sample has different column order or drops a variable, the result is reordered to match the original. If states are missing, the iteration is marked as failed.

**Key**: Multilevel decomposition is applied to each bootstrap sample, not just to the original. This correctly captures the uncertainty in the decomposition.

---

## 7. Statistics: `.compute_bootstrap_stats()` (lines 410-469)

Input: `boot_matrices` (iter x nbins), `original_matrix`, `states`, `directed`, etc.

### Mean and SD (lines 417-418)

```r
weights_mean <- colMeans(boot_matrices, na.rm = TRUE)
weights_sd   <- apply(boot_matrices, 2, sd, na.rm = TRUE)
```

### Percentile CIs (lines 421-424)

```r
ci_lower <- apply(boot_matrices, 2, quantile, probs = ci_level / 2, na.rm = TRUE)
ci_upper <- apply(boot_matrices, 2, quantile, probs = 1 - ci_level / 2, na.rm = TRUE)
```

### P-values (lines 427-446)

```r
valid_rows <- rowSums(is.na(boot_matrices)) == 0
bm <- boot_matrices[valid_rows, , drop = FALSE]
n_valid <- nrow(bm)
```

**Stability inference** (lines 433-441):

```r
cr_low  <- pmin(orig_flat * consistency_range[1], orig_flat * consistency_range[2])
cr_high <- pmax(orig_flat * consistency_range[1], orig_flat * consistency_range[2])
below   <- sweep(bm, 2, cr_low, "<")
above   <- sweep(bm, 2, cr_high, ">")
p_counts <- colSums(below | above)
p_values <- (p_counts + 1) / (n_valid + 1)
```

**Key**: `pmin`/`pmax` handles negative weights correctly (e.g., if `orig = -0.3`, then `cr_low = -0.3 * 1.25 = -0.375` and `cr_high = -0.3 * 0.75 = -0.225`).

**Key**: `(p_counts + 1) / (n_valid + 1)` is the conservative p-value formula (adds pseudocount to avoid p = 0).

**Threshold inference** (lines 443-446):

```r
below_thresh <- sweep(abs(bm), 2, edge_threshold, "<")
p_counts <- colSums(below_thresh)
p_values <- (p_counts + 1) / (n_valid + 1)
```

### Significant mask (lines 453-455)

```r
p_mat    <- to_mat(p_values)
sig_mask <- (p_mat < ci_level) * 1
sig_mat  <- original_matrix * sig_mask
```

Multiplying by the mask zeros out non-significant edges while preserving original weights for significant ones.

---

## 8. Summary: `.build_bootstrap_summary()` (lines 474-501)

Builds a long-format data frame using `data.table::data.table()`.

```r
dt <- data.table::data.table(
  from     = rep(states, each = n),
  to       = rep(states, times = n),
  weight   = as.vector(t(original_matrix)),
  mean     = as.vector(t(stats$mean)),
  sd       = as.vector(t(stats$sd)),
  p_value  = as.vector(t(stats$p_values)),
  sig      = as.vector(t(stats$p_values)) < ci_level,
  ci_lower = as.vector(t(stats$ci_lower)),
  ci_upper = as.vector(t(stats$ci_upper))
)
```

**CRITICAL**: `as.vector(t(matrix))` — the `t()` transpose is essential for row-major flattening. `as.vector()` on a matrix in R flattens column-major. Transposing first gives row-major order, which matches `from = rep(states, each = n)` (each source state paired with all targets).

**Filtering** (lines 494-498):

```r
if (directed) {
  dt <- dt[weight != 0 & from != to]   # non-zero, non-self-loop
} else {
  dt <- dt[weight != 0 & from < to]    # upper triangle only
}
```

---

## 9. Pruned Model (lines 179-196)

```r
pruned_matrix <- stats$significant
pruned_edges  <- .extract_edges_from_matrix(pruned_matrix, directed = directed)

model <- list(
  matrix = pruned_matrix,
  nodes = states,
  directed = directed,
  method = method,
  params = params,
  scaling = scaling,
  threshold = threshold,
  n_nodes = n_states,
  n_edges = nrow(pruned_edges),
  edges = pruned_edges,
  level = level
)
class(model) <- "netobject"
```

The pruned model is a `netobject` with non-significant edges zeroed out. It can be used directly for downstream analysis.

---

## 10. Return Value (lines 198-219)

```r
result <- list(
  original          = original,          # netobject
  mean              = stats$mean,        # matrix
  sd                = stats$sd,          # matrix
  p_values          = stats$p_values,    # matrix
  significant       = stats$significant, # matrix
  ci_lower          = stats$ci_lower,    # matrix
  ci_upper          = stats$ci_upper,    # matrix
  cr_lower          = stats$cr_lower,    # matrix (stability) or NULL
  cr_upper          = stats$cr_upper,    # matrix (stability) or NULL
  summary           = summary_df,        # data.frame
  model             = model,             # netobject (pruned)
  method            = method,
  params            = params,
  iter              = iter,
  ci_level          = ci_level,
  inference         = inference,
  consistency_range = consistency_range,
  edge_threshold    = edge_threshold
)
class(result) <- "saqr_bootstrap"
```

---

## 11. Auto Edge Threshold (lines 145-153)

When `inference = "threshold"` and `edge_threshold` is NULL, it's auto-set to the 10th percentile of absolute non-zero original weights:

```r
if (inference == "threshold" && is.null(edge_threshold)) {
  abs_weights <- abs(as.vector(original$matrix))
  nz_weights <- abs_weights[abs_weights > 0]
  edge_threshold <- if (length(nz_weights) > 0) {
    quantile(nz_weights, probs = 0.10)
  } else {
    0
  }
}
```

---

## 12. Performance Characteristics

| Aspect | Transition fast path | Association path |
|---|---|---|
| Per-iteration cost | `colSums()` (C-level) + inline post-process | Full estimator call |
| Memory | `n_seq x nbins` precompute matrix | No persistent storage |
| Speed vs tna::bootstrap() | ~2.8x faster (2000 sequences, 500 iter) | N/A (different methods) |
| Why faster | `colSums` on 2D matrix vs `apply` on 3D array | — |
| Format requirement | Wide only (errors on long) | Any (estimator handles) |
